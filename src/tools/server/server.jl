include("../../chunkedgraphs/ChunkedGraphs.jl")
using ChunkedGraphs

import CloudVolume
import MbedTLS
import HttpServer
import HttpCommon

import Memento
import JSON

import Base.run

struct ServerConfig
	server::HttpServer.Server
	host::IPAddr
	port::Integer
	ssl::Tuple{MbedTLS.CRT, MbedTLS.PKContext}
end

mutable struct ChunkedGraphServer
	cgraph::ChunkedGraph{Chunk}
	logger::Memento.Logger
	serverconfig::ServerConfig
	watershed::CloudVolume.CloudVolumeWrapper

	rg2cg::Dict{ChunkID, Dict{UInt64, Label}}
	cg2rg::Dict{Label, UInt64}

	function ChunkedGraphServer(config::AbstractString)
		cgs = new()

		# Used when interacting with original watershed layer
		cgs.rg2cg = Dict{ChunkID, Dict{UInt64, Label}}()
		cgs.cg2rg = Dict{Label, UInt64}()

		# Read config file
		isfile(config) || error("Chunked Graph Server configuration not found: $config")
		settings = JSON.parsefile(config)
		basedir = dirname(config)

		# Prepare dir structure
		mkpath(abspath(settings["graphpath"], basedir))
		mkpath(abspath(settings["logpath"], basedir))
		mkpath(abspath(settings["certpath"], basedir))

		# Load/Initialize Chunked Graph
		cgs.cgraph = ChunkedGraph(abspath(settings["graphpath"], basedir))
		# gc_enable(false)
		# @time for f in filter(s->ismatch(r".*\.chunk", s), readdir(abspath(settings["graphpath"], basedir)))
		# 	m = match(r"(\d+)_(\d+)_(\d+)_(\d+)\..*", f)
		# 	id = tochunkid(map(x->parse(UInt32, x), m.captures)...)
		# 	if tolevel(id) >= 3
		# 		getchunk!(cgs.cgraph, id)
		# 	end
		# end
		# gc_enable(true)

		# CloudVolume to fetch watershed layer if necessary
		cgs.watershed = CloudVolume.CloudVolumeWrapper(settings["cloudpath"]; bounded=false, fill_missing=true)

		# Set up logger
		cgs.logger = Memento.getlogger("Graph Operations")
		Memento.setlevel!(cgs.logger, "info")
		push!(cgs.logger, Memento.DefaultHandler(
				joinpath(abspath(settings["logpath"], basedir), "graph.log"),
				Memento.DefaultFormatter("[{date} | {level} | {name}]: {msg}")))

		# Generate a certificate and key if they do not exist
		if !isfile(joinpath(abspath(settings["certpath"], basedir), "server.crt"))
			@static if is_unix()
				run(`openssl req -x509 -nodes -days 365 -newkey rsa:2048 -keyout
					$(joinpath(abspath(settings["certpath"], basedir), "server.key")) -out $(joinpath(abspath(settings["certpath"], basedir), "server.crt"))`)
			end
		end
		cert = MbedTLS.crt_parse_file(joinpath(abspath(settings["certpath"], basedir), "server.crt"))
		key = MbedTLS.parse_keyfile(joinpath(abspath(settings["certpath"], basedir), "server.key"))

		# Prepare endpoints
		headers = HttpCommon.headers()
		headers["Access-Control-Allow-Origin"]= "*"
		headers["Access-Control-Allow-Headers"]= "Origin, X-Requested-With, Content-Type, Accept"
		headers["Access-Control-Allow-Methods"]= "POST, GET, OPTIONS"

		http = HttpServer.HttpHandler() do req::HttpServer.Request, res::HttpServer.Response
			if req.method == "OPTIONS"
				return HttpServer.Response(UInt8[], headers)
			elseif ismatch(r"/1.0/segment/(\d+)/children/?", req.resource) && req.method == "GET"
				return HttpServer.Response(handle_children(cgs, match(r"/1.0/segment/(\d+)/children/?", req.resource).captures[1]), headers)
			elseif ismatch(r"/1.0/segment/(\d+)/leaves/?", req.resource) && req.method == "GET"
				return HttpServer.Response(handle_leaves(cgs, match(r"/1.0/segment/(\d+)/leaves/?(?:\?(.*))?", req.resource).captures...), headers)
			elseif ismatch(r"/1.0/graph/root/?", req.resource) && req.method == "POST"
				return HttpServer.Response(handle_root(cgs, req.data), headers)
			elseif ismatch(r"/1.0/graph/merge/?", req.resource) && req.method == "POST"
				return HttpServer.Response(handle_merge(cgs, req.data), headers)
			elseif ismatch(r"/1.0/graph/split/?", req.resource) && req.method == "POST"
				return HttpServer.Response(handle_split(cgs, req.data), headers)
			elseif ismatch(r"/1.0/graph/save/?", req.resource) && req.method == "POST"
				return HttpServer.Response(handle_save(cgs), headers)
			else
				println("could not parse $(req.resource)")
				return HttpServer.Response(400)
			end
		end
		http.events["listen"] = (saddr) -> println("Running on https://$saddr (Press CTRL+C to quit)")

		server = HttpServer.Server(http)
		host = getaddrinfo(settings["host"])
		port = settings["port"]
		ssl = (cert, key)

		cgs.serverconfig = ServerConfig(server, host, port, ssl)
		return cgs
	end
end

function run(cgs::ChunkedGraphServer)
	run(cgs.serverconfig.server; host = cgs.serverconfig.host, port = cgs.serverconfig.port, ssl = cgs.serverconfig.ssl)
end

function simple_print(x::Array)
	string('[', map(n->"$(n),", x)..., ']')
end

# TODO: Eventually we'll have to clear old labels
function loadlabels!(cgs::ChunkedGraphServer, chunkid::ChunkID)
	if !haskey(cgs.rg2cg, chunkid)
		cg2rg = Dict{Label, UInt64}(reinterpret(Pair{Label, UInt64}, read("$(cgs.cgraph.path)cg2rg_$(ChunkedGraphs.stringify(chunkid)).bin")))
		merge!(cgs.cg2rg, cg2rg)
		cgs.rg2cg[chunkid] = map(reverse, cg2rg)
	end
end

function pos2cg(cgs::ChunkedGraphServer, pos::Tuple{Integer, Integer, Integer},
		cgidhint::Label, default::Union{Label,Void} = nothing)
	cgidroot = root!(cgs.cgraph, getvertex!(cgs.cgraph, cgidhint))

	padding = (16, 16, 2)
	voxelres = cgs.watershed.val[:scale]["resolution"]
	cutoutrange = @. range.(pos - padding, 2 * padding)

	svids = cgs.watershed[cutoutrange...]
	worldpositions = collect(CartesianRange(cutoutrange))

	const sortkey = sortperm(collect(Base.Iterators.flatten(
		map(i -> sum(@. (voxelres * (i.I - padding - 1)) ^ 2), CartesianRange(size(svids)))
	)))

	checked = Set{UInt64}()
	for i in sortkey
		if svids[i] == 0 || svids[i] in checked
			continue
		end

		cgid = rg2cg(cgs, svids[i], worldpositions[i].I)

		if root!(cgs.cgraph, getvertex!(cgs.cgraph, cgid)) == cgidroot
			return cgid
		end

		push!(checked, svids[i])
	end
	return default
end

function rg2cg(cgs::ChunkedGraphServer, rgid::UInt64, pos::Tuple{Integer,Integer,Integer},
		default::Union{Label,Void} = nothing)
	chunkid = ChunkedGraphs.world_to_chunk(pos...)
	loadlabels!(cgs, chunkid)
	return get(cgs.rg2cg[chunkid], rgid, default)
end

function cg2rg(cgs::ChunkedGraphServer, cgid::Label,
		default::Union{Label,Void} = nothing)
	loadlabels!(cgs, tochunkid(cgid))
	return get(cgs.cg2rg, cgid, default)
end

# Input IDs should be a ChunkedGraph ID
# Output IDs are Neuroglancer supervoxel IDs
function handle_leaves(cgs::ChunkedGraphServer, cgid::AbstractString, query::Union{AbstractString, Void})
	Memento.debug(cgs.logger, "handle_leaves($cgid)")

	cgid = parse(Label, cgid)
	bbox = nothing

	if query !== nothing
		matches = match(r"bounds=(\d+)-(\d+)_(\d+)-(\d+)_(\d+)-(\d+)", query)
		if matches !== nothing
			bounds = map(x -> parse(Int, x), matches.captures)
			chunk_min = fld(bounds[1], ChunkedGraphs.CHUNK_SIZE[1]),
						fld(bounds[3], ChunkedGraphs.CHUNK_SIZE[2]),
						fld(bounds[5], ChunkedGraphs.CHUNK_SIZE[3])
			chunk_max = fld(bounds[2] - 1, ChunkedGraphs.CHUNK_SIZE[1]),
						fld(bounds[4] - 1, ChunkedGraphs.CHUNK_SIZE[2]),
						fld(bounds[6] - 1, ChunkedGraphs.CHUNK_SIZE[3])
			bbox = (chunk_min[1]:chunk_max[1], chunk_min[2]:chunk_max[2], chunk_min[3]:chunk_max[3])
		end
	end

	ancestor = getvertex!(cgs.cgraph, cgid)
	segments = leaves!(cgs.cgraph, ancestor, 1, bbox)

	println("$(now()): selected $(length(segments)) segments with ancestor $(ancestor.label) in region $(bbox)")
	s = collect(Set{UInt64}(cg2rg(cgs, x) for x in segments))

	return reinterpret(UInt8, s)
end

# Input IDs should be a ChunkedGraph ID
# Output IDs are ChunkedGraph IDs
function handle_children(cgs::ChunkedGraphServer, cgid::AbstractString)
	Memento.debug(cgs.logger, "handle_children($cgid)")

	cgid = parse(UInt64, cgid)

	v = getvertex!(cgs.cgraph, cgid)

	if tolevel(v) == 1
		# Lvl 1, convert ChunkedGraph Supervoxel to Neuroglancer supervoxel
		s = UInt64[cg2rg(cgs, v.label)]
		println("$(now()): handle_children - v: $(v.label), (Level $(tolevel(v))), - RG equivalent: $(cg2rg(cgs, v.label))")
	elseif tolevel(v) == 2
		# Lvl 2, convert child nodes (ChunkedGraph supervoxel) to Neuroglancer supervoxel
		s = map(x->cg2rg(cgs, x), v.children)
		println("$(now()): handle_children - v: $(v.label), (Level $(tolevel(v))), - children: $(simple_print(map(x->cg2rg(cgs, x), v.children)))")
	else
		# Lvl 3+
		# s = UInt64[child for child in v.children]
		s = UInt64[child for child in leaves!(cgs.cgraph, v, 2)] # J's hack to skip the middle layers and jump right to the pre-meshed lower level agglomeration.
		println("$(now()): handle_children - v: $(v.label), (Level $(tolevel(v))), - children: $(simple_print([child for child in v.children]))")
	end

	return reinterpret(UInt8, s)
end

# Input IDs should be a Neuroglancer supervoxel ID
# Output ID is a ChunkedGraph ID
function handle_root(cgs::ChunkedGraphServer, data::Vector{UInt8})
	# FIXME: What if this is called with a valid CG ID?
	Memento.debug(cgs.logger, "handle_root($(String(data)))")

	parsed = JSON.parse(String(data))
	id = parse(UInt64, parsed[1])
	pos = (parsed[2], parsed[3], parsed[4])
	print("$(now()): Root for segment $(id) at position $(pos): ")

	id = rg2cg(cgs, id, pos)
	if id isa Void
		println("ID not found in chunk - wrong position?")
		return UInt8[]
	end
	root_vertex = root!(cgs.cgraph, getvertex!(cgs.cgraph, id))

	println("$(root_vertex.label)")
	return reinterpret(UInt8, [root_vertex.label])
end

# Input IDs should be ChunkedGraph IDs
# Output IDs are ChunkedGraph IDs
function handle_split(cgs::ChunkedGraphServer, data::Vector{UInt8})
	Memento.info(cgs.logger, "handle_split($(String(data)))")

	parsed = JSON.parse(String(data))
	parsedsources = map(x -> (parse(UInt64, x[1]), (x[2], x[3], x[4])), parsed["sources"])
	parsedsinks   = map(x -> (parse(UInt64, x[1]), (x[2], x[3], x[4])), parsed["sinks"])

	sources = map(x -> pos2cg(cgs, x[2], rg2cg(cgs, x[1], x[2], x[1])), parsedsources)
	sources = unique(convert(Vector{UInt64}, filter(x -> !isa(x, Void), sources)))

	sinks   = map(x -> pos2cg(cgs, x[2], rg2cg(cgs, x[1], x[2], x[1])), parsedsinks)
	sinks   = unique(convert(Vector{UInt64}, filter(x -> !isa(x, Void), sinks)))

	if isempty(sources) || isempty(sinks)
		println("Empty source or sink.")
		return UInt8[]
	end

	if !isempty(intersect(sources, sinks))
		println("Source and sink are the same")
		return UInt8[]
	end

	cuts = mincut!(cgs.cgraph, sources, sinks)

	for e in cuts
		delete_atomic_edge!(cgs.cgraph, e)
	end

	update!(cgs.cgraph)

	root_labels = Set{UInt64}()
	for e in cuts
		push!(root_labels, root!(cgs.cgraph, getvertex!(cgs.cgraph, e.u)).label)
		push!(root_labels, root!(cgs.cgraph, getvertex!(cgs.cgraph, e.v)).label)
	end

	root_labels = collect(root_labels) #Vector{UInt64}(map(cgid -> tolevel(cgid) == 1 ? cg2rg(cgs, cgid) : cgid, collect(root_labels)))

	println("$(now()): Split $(sources) and $(sinks) => $(simple_print(root_labels))")
	return reinterpret(UInt8, root_labels)
end

# Input IDs should be ChunkedGraph IDs
# Output ID is a ChunkedGraph ID
function handle_merge(cgs::ChunkedGraphServer, data::Vector{UInt8})
	Memento.info(cgs.logger, "handle_merge($(String(data)))")

	parsed = JSON.parse(String(data))
	parsedsegments = map(x -> (parse(UInt64, x[1]), (x[2], x[3], x[4])), parsed)
	
	segments = map(x -> pos2cg(cgs, x[2], rg2cg(cgs, x[1], x[2], x[1])), parsedsegments)
	segments = unique(convert(Vector{UInt64}, filter(x -> !isa(x, Void), segments)))

	@assert length(segments) == 2

	add_atomic_edge!(cgs.cgraph, AtomicEdge(segments[1], segments[2]))

	update!(cgs.cgraph)

	root = root!(cgs.cgraph, getvertex!(cgs.cgraph, segments[1]))

	println("$(now()): Merged $(cg2rg(cgs, segments[1])) and $(cg2rg(cgs, segments[2])) => $(root.label)")
	return reinterpret(UInt8, [root.label])
end

function handle_save(cgs::ChunkedGraphServer)
	Memento.debug(cgs.logger, "handle_save()")
	update!(cgs.cgraph)
	save!(cgs.cgraph)
	
	return UInt8[]
end

cgs = ChunkedGraphServer("server.conf")
#run(cgs)
