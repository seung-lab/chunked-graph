include("../../chunkedgraphs/ChunkedGraphs.jl")
using ChunkedGraphs
using PyCall
import JSON

@pyimport taskqueue
@pyimport igneous.tasks as tasks

@inline function stringify(x::UnitRange{T}) where T <: Integer
	return "$(x.start)-$(x.stop)"
end

@inline function stringify(x::Tuple{UnitRange{T}, UnitRange{T}, UnitRange{T}}) where T <: Integer
	return join(map(stringify, x), '_')
end

mutable struct ChunkedGraphMesher
	cgraph::ChunkedGraph
	rg2cg::Dict{ChunkID, Dict{UInt64, Label}}
	cg2rg::Dict{Label, UInt64}
	settings::Dict{String, Any}
	cloudpath::AbstractString
	mip::Int
	simplification_factor::Real
	max_simplification_error::Real
	low_padding::Int
	high_padding::Int

	function ChunkedGraphMesher(config::AbstractString)
		cgm = new()

		cgm.rg2cg = Dict{ChunkID, Dict{UInt64, Label}}()
		cgm.cg2rg = Dict{Label, UInt64}()

		# Read config file
		isfile(config) || error("Chunked Graph Mesher configuration not found: $config")
		settings = JSON.parsefile(config)
		basedir = dirname(config)

		# Load/Initialize Chunked Graph
		cgm.cgraph = ChunkedGraph(abspath(settings["graphpath"], basedir))

		# Settings
		cgm.cloudpath = settings["cloudpath"]
		cgm.mip = settings["mip"]
		cgm.simplification_factor = settings["simplification_factor"]
		cgm.max_simplification_error = settings["max_simplification_error"]
		cgm.low_padding = settings["low_padding"]
		cgm.high_padding = settings["high_padding"]

		return cgm
	end
end

# TODO: Eventually we'll have to clear old labels
function loadlabels!(cgm::ChunkedGraphMesher, chunkid::ChunkID)
	if !haskey(cgm.rg2cg, chunkid)
		cg2rg = Dict{Label, UInt64}(reinterpret(Pair{Label, UInt64}, read("$(cgm.cgraph.path)cg2rg_$(ChunkedGraphs.stringify(chunkid)).bin")))
		merge!(cgm.cg2rg, cg2rg)
		cgm.rg2cg[chunkid] = map(reverse, cg2rg)
	end
end

function cg2rg(cgm::ChunkedGraphMesher, cgid::Label,
		default::Union{Label,Void} = nothing)
	loadlabels!(cgm, tochunkid(cgid))
	return get(cgm.cg2rg, cgid, default)
end

function rg2cg(cgm::ChunkedGraphMesher, rgid::UInt64, pos::Tuple{Integer,Integer,Integer},
		default::Union{Label,Void} = nothing)
	chunkid = ChunkedGraphs.world_to_chunkid(cgm.cgraph, pos)
	loadlabels!(cgm, chunkid)
	return get(cgm.rg2cg[chunkid], rgid, default)
end

function mesh!(cgm::ChunkedGraphMesher, c::Chunk)
	level = tolevel(c)
	if level == 1
		tic();
		print("Meshing Chunk $(num2hex(c.id)): ")
		vertexlist = Set{Vertex}()
		# TODO: Ideally we should include vertices in the overlapping region
		for v in values(c.vertices)
			if v.parent !== NULL_LABEL
				push!(vertexlist, getvertex!(c.cgraph, v.parent))
			end
		end

		remap = Dict{Label, Label}()
		for v in vertexlist
			for l in leaves!(c.cgraph, v)
				remap[cg2rg(cgm, l)] = v.label
			end
		end
		println("$(length(remap)) supervoxel, $(length(vertexlist)) agglomerations ... ")

		task = tasks.MeshTask(
			div.(c.cgraph.CHUNKSIZE, (8, 8, 1)),
			div.(start.(ChunkedGraphs.chunkid_to_ranges(c.cgraph, c.id)), (8, 8, 1)),
			cgm.cloudpath,
			lod = 0,
			mip = cgm.mip,
			simplification_factor = cgm.simplification_factor,
			max_simplification_error = cgm.max_simplification_error,
			remap_table = remap,
			generate_manifests = true,
			low_padding = cgm.low_padding,
			high_padding = cgm.high_padding
		)
		task[:execute]()
		println("Done ($(toq()) s). Check $(first(vertexlist).label), for example.")
	end
end

function mesh!(cgm::ChunkedGraphMesher)
	@time for f in filter(s->ismatch(r".*\.chunk", s), readdir(cgm.cgraph.path))
		m = match(r"(\d+)_(\d+)_(\d+)_(\d+)\..*", f)
		id = tochunkid(map(x->parse(UInt32, x), m.captures)...)
		if tolevel(id) == 1
			mesh!(cgm, getchunk!(cgm.cgraph, id))
			return
		end
	end

end

cgm = ChunkedGraphMesher("./pinky40.conf")
