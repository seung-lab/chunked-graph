#TODO: Cleanup
push!(LOAD_PATH, dirname(@__FILE__))
include("./chunkedgraphs/ChunkedGraphs.jl")
using ChunkedGraphs

import MbedTLS
import HttpServer
import HttpCommon

#import Logging # TODO: Broken in 0.6
import JSON

rel(p::String) = joinpath(dirname(@__FILE__), p)

# Prepare dir structure
settings = JSON.parsefile(rel("server.conf"))
@static if is_unix()
	run(`mkdir -p $(rel(settings["graphpath"]))`)
	run(`mkdir -p $(rel(settings["logpath"]))`)
	run(`mkdir -p $(rel(settings["certpath"]))`)
end

# Basic logging
#@Logging.configure(level = DEBUG)
#Logging.configure(filename = joinpath(rel(settings["logpath"]), "graph.log"))

# Generate a certificate and key if they do not exist
if !isfile(joinpath(rel(settings["certpath"]), "server.crt"))
	@static if is_unix()
		run(`openssl req -x509 -nodes -days 365 -newkey rsa:2048 -keyout
			$(joinpath(rel(settings["certpath"]), "server.key")) -out $(joinpath(rel(settings["certpath"]), "server.crt"))`)
	end
end

# Load top layers of Chunked Graph

cgraph = ChunkedGraph(rel(settings["graphpath"]), settings["cloudpath"])
gc_enable(false)
@time for f in filter(s->ismatch(r".*\.chunk", s), readdir(expanduser(rel(settings["graphpath"]))))
	m = match(r"(\d+)_(\d+)_(\d+)_(\d+)\..*", f)
	id = tochunkid(map(x->parse(UInt32, x), m.captures)...)
	if tolevel(id) >= 3
		getchunk!(cgraph, id)
	end
end
gc_enable(true)

rg2cg = Dict{UInt64, Label}(reinterpret(Pair{UInt64, Label}, read(open(joinpath(rel(settings["graphpath"]), "rg_to_cg.bin"), "r"))))
cg2rg = Dict{Label, UInt64}(reinterpret(Pair{Label, UInt64}, read(open(joinpath(rel(settings["graphpath"]), "cg_to_rg.bin"), "r"))))

function simple_print(x::Array)
	string('[', map(n->"$(n),", x)..., ']')
end

function handle_leaves(id::AbstractString, query::Union{AbstractString, Void})
	#@Logging.debug("handle_leaves($id)")
	id = parse(UInt64, id)
	if tochunkid(id) == 0 # Lvl 1, a neuroglancer supervoxel, need to lookup chunk id
		id = rg2cg[id]
	end

	bbox = (0:typemax(Int), 0:typemax(Int), 0:typemax(Int))

	if query !== nothing
		matches = match(r"bounds=(\d+)-(\d+)_(\d+)-(\d+)_(\d+)-(\d+)", query)
		if matches !== nothing
			bounds = map(x->parse(Int, x), matches.captures)
			chunk_min = fld(bounds[1], ChunkedGraphs.CHUNK_SIZE[1]),
						fld(bounds[3], ChunkedGraphs.CHUNK_SIZE[2]),
						fld(bounds[5], ChunkedGraphs.CHUNK_SIZE[3])
			chunk_max = fld(bounds[2] - 1, ChunkedGraphs.CHUNK_SIZE[1]),
						fld(bounds[4] - 1, ChunkedGraphs.CHUNK_SIZE[2]),
						fld(bounds[6] - 1, ChunkedGraphs.CHUNK_SIZE[3])
			bbox = (chunk_min[1]:chunk_max[1], chunk_min[2]:chunk_max[2], chunk_min[3]:chunk_max[3])
		end
	end

	ancestor = getvertex!(cgraph, id)
	segments = leaves!(cgraph, ancestor, 1, bbox)

	println("$(now()): selected $(length(segments)) segments with ancestor $(ancestor.label) in region $(bbox)")
	s = collect(Set{UInt64}(cg2rg[x] for x in segments))

	return HttpServer.Response(reinterpret(UInt8, s), headers)
end

function handle_root(id::AbstractString)
	#@Logging.debug("handle_root($id)")
	id = parse(UInt64, id)
	print("$(now()): Root for segment $(id): ")

	if tochunkid(id) == 0 # Lvl 1, a neuroglancer supervoxel, need to lookup chunk id
		id = rg2cg[id]
	end

	root_vertex = [root!(cgraph, getvertex!(cgraph, id))][1]
	println("$(root_vertex.label)")

	return HttpServer.Response(reinterpret(UInt8,[root_vertex.label]),headers)
end

function handle_children(id::AbstractString)
	#@Logging.debug("handle_children($id)")
	id = parse(UInt64, id)

	if tochunkid(id) == 0 # Lvl 1, a neuroglancer supervoxel, need to lookup chunk id
		id = rg2cg[id]
	end

	v = getvertex!(cgraph, id)

	if tolevel(v) == 1 # Lvl 1, a neuroglancer supervoxel, no children
		s = UInt64[]
		println("$(now()): handle_children - v: $(v.label), (Level $(tolevel(v)))")
	elseif tolevel(v) == 2 # Lvl 2, children are neuroglancer supervoxel, need to trim the chunk ids
		s = map(x->cg2rg[x], v.children)
		println("$(now()): handle_children - v: $(v.label), (Level $(tolevel(v))), - children: $(simple_print(map(x->cg2rg[x], v.children)))")
	else
		#s = UInt64[child for child in v.children]
		s = UInt64[child for child in leaves!(cgraph,v,2)] # J's hack to skip the middle layers and jump right to the pre-meshed lower level agglomeration.
		println("$(now()): handle_children - v: $(v.label), (Level $(tolevel(v))), - children: $(simple_print([child for child in v.children]))")
	end

	return HttpServer.Response(reinterpret(UInt8,s),headers)
end

function handle_split(data::Vector{UInt8})
	#@Logging.info("handle_split($(String(data)))")
	parsed = JSON.parse(String(data))

	sources = unique(convert(Vector{UInt64}, filter(y->y != nothing, map(x->getsupervoxelat(cgraph, parse(UInt64, x[1]), (x[2], x[3], x[4])), parsed["sources"]))))
	sinks = unique(convert(Vector{UInt64}, filter(y->y != nothing, map(x->getsupervoxelat(cgraph, parse(UInt64, x[1]), (x[2], x[3], x[4])), parsed["sinks"]))))

	if isempty(sources) || isempty(sinks)
		println("Empty source or sink.")
		return HttpServer.Response(UInt8[], headers)
	end
	if !isempty(intersect(sources, sinks))
		println("Source and sink are the same")
		return HttpServer.Response(UInt8[], headers)
	end

	cuts = mincut!(cgraph, sources, sinks)
	for e in cuts
		delete_atomic_edge!(cgraph,e)
	end
	update!(cgraph)

	root_labels = Set{UInt64}()
	for e in cuts
		push!(root_labels, root!(cgraph, getvertex!(cgraph, e.u)).label)
		push!(root_labels, root!(cgraph, getvertex!(cgraph, e.v)).label)
	end

	root_labels = Array{UInt64}(map(x->tolevel(tochunkid(x)) == 1 ? cg2rg[x] : x, collect(root_labels)))

	println("$(now()): Split $(sources) and $(sinks) => $(simple_print(root_labels))")
	return HttpServer.Response(reinterpret(UInt8, root_labels), headers)
end

function handle_merge(data::Vector{UInt8})
	#@Logging.info("handle_merge($(String(data)))")
	parsed = JSON.parse(String(data))

	segments = unique(convert(Vector{UInt64}, filter(y->y != nothing, map(x->getsupervoxelat(cgraph, parse(UInt64, x[1]), (x[2], x[3], x[4])), parsed))))
	@assert length(segments) == 2

	add_atomic_edge!(cgraph, AtomicEdge(segments[1], segments[2]))
	update!(cgraph)

	root = root!(cgraph, getvertex!(cgraph, segments[1]))
	println("$(now()): Merged $(cg2rg[segments[1]]) and $(cg2rg[segments[2]]) => $(root.label)")
	
	return HttpServer.Response(reinterpret(UInt8, [root.label]), headers)
end

function handle_save()
	#@Logging.info("handle_save()")
	update!(cgraph)
	save!(cgraph)
	return HttpServer.Response(UInt8[], headers)
end

headers = HttpCommon.headers()

headers["Access-Control-Allow-Origin"]= "*"
headers["Access-Control-Allow-Headers"]= "Origin, X-Requested-With, Content-Type, Accept"
headers["Access-Control-Allow-Methods"]= "POST, GET, OPTIONS"

http = HttpServer.HttpHandler() do req::HttpServer.Request, res::HttpServer.Response
	if req.method == "OPTIONS"
		return HttpServer.Response(UInt8[], headers)
	elseif ismatch(r"/1.0/segment/(\d+)/root/?", req.resource) && req.method == "GET"
		return handle_root(match(r"/1.0/segment/(\d+)/root/?", req.resource).captures[1])
	elseif ismatch(r"/1.0/segment/(\d+)/children/?", req.resource) && req.method == "GET"
		return handle_children(match(r"/1.0/segment/(\d+)/children/?", req.resource).captures[1])
	elseif ismatch(r"/1.0/segment/(\d+)/leaves/?", req.resource) && req.method == "GET"
		return handle_leaves(match(r"/1.0/segment/(\d+)/leaves/?(?:\?(.*))?", req.resource).captures...)
	elseif ismatch(r"/1.0/graph/merge/?", req.resource) && req.method == "POST"
		return handle_merge(req.data)
	elseif ismatch(r"/1.0/graph/split/?", req.resource) && req.method == "POST"
		return handle_split(req.data)
	elseif ismatch(r"/1.0/graph/save/?", req.resource) && req.method == "POST"
		return handle_save()
	else
		println("could not parse $(req.resource)")
		return HttpServer.Response(400)
	end
end
http.events["listen"] = (saddr) -> println("Running on https://$saddr (Press CTRL+C to quit)")

server = HttpServer.Server(http)

cert = MbedTLS.crt_parse_file(joinpath(rel(settings["certpath"]), "server.crt"))
key = MbedTLS.parse_keyfile(joinpath(rel(settings["certpath"]), "server.key"))

#run(server, host=getaddrinfo(settings["host"]), port=settings["port"], ssl=(cert, key))
