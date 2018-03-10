push!(LOAD_PATH, dirname(@__FILE__))
include("../chunkedgraphs/ChunkedGraphs.jl")
using ChunkedGraphs
using CloudVolume
using CodecZstd
using DataFrames
using Memoize
using MySQL
using OffsetArrays

import PyCall

const RG_CHUNKSIZE = (512, 512, 512)
const CG_CHUNKSIZE = (512, 512, 64)
const mysql_conn = MySQL.connect("127.0.0.1", "root", readline("/secrets/mysql"); db="relabeling")

mutable struct RegionGraphWrapper
	storage::StorageWrapper
	chunksize::Tuple{Unsigned, Unsigned, Unsigned}
	offset::Tuple{Integer, Integer, Integer}
	maxlevel::Unsigned
end

mutable struct RelabeledChunk{T <: Integer}
	watershed_original::OffsetArray{T, 3}
	watershed_relabeled::OffsetArray{Label, 3}
	rg_to_cg_complete::Dict{ChunkID, Dict{T, Label}}
	rg_to_cg_boundary::Dict{ChunkID, Dict{T, Label}}
end

function RelabeledChunk(watershed::OffsetArray{T, 3}) where T <: Integer
	RelabeledChunk{T}(watershed,
				   OffsetArray(Label, indices(watershed)...),
				   Dict{ChunkID, Dict{T, Label}}(),
				   Dict{ChunkID, Dict{T, Label}}()
	)
end

@inline function drop_last(x::UnitRange{T}) where T <: Integer
	return x.start : x.stop - 1
end

@inline function stringify(x::UnitRange{T}) where T <: Integer
	return "$(x.start)-$(x.stop)"
end

@inline function stringify(x::Tuple{UnitRange{T}, UnitRange{T}, UnitRange{T}}) where T <: Integer
	return join(map(stringify, x), '_')
end

@inline function toslice(s::AbstractString)
	m = match(r"(\d+)-(\d+)_(\d+)-(\d+)_(\d+)-(\d+)", s)
	bounds = map(x -> parse(Int, x), m.captures)
	return (bounds[1] : bounds[2], bounds[3] : bounds[4], bounds[5] : bounds[6])
end

function load_labels(watershed::OffsetArray{T, 3}) where T <: Integer
	voxel_range = indices(watershed)[1:3]
	chunk_range = map(x -> div.(x, CG_CHUNKSIZE), (minimum.(voxel_range), maximum.(voxel_range)))
	res = MySQL.query(mysql_conn, """
		SELECT c.x, c.y, c.z, m.regiongraph_sv, m.chunkedgraph_sv
		FROM chunks c JOIN mappings m ON c.id = m.chunk_id
		WHERE c.x BETWEEN $(chunk_range[1][1]) AND $(chunk_range[2][1] + 1) AND
		      c.y BETWEEN $(chunk_range[1][2]) AND $(chunk_range[2][2] + 1) AND
		      c.z BETWEEN $(chunk_range[1][3]) AND $(chunk_range[2][3] + 1);
	""", DataFrame)

	chunk_labels = Dict{ChunkID, Dict{T, Label}}()
	for row in eachrow(res)
		chunkid = tochunkid(1, row[:x], row[:y], row[:z])
		if !haskey(chunk_labels, chunkid)
			chunk_labels[chunkid] = Dict{T, Label}(row[:regiongraph_sv] => row[:chunkedgraph_sv])
		else
			push!(chunk_labels[chunkid], row[:regiongraph_sv] => row[:chunkedgraph_sv])
		end
	end
	return chunk_labels::Dict{ChunkID, Dict{T, Label}}
end

function save_labels(chunk_labels::Dict{ChunkID, Dict{T, Label}}) where T <: Integer
	if length(chunk_labels) > 0
		MySQL.execute!(mysql_conn,
			"INSERT IGNORE INTO chunks (x,y,z) VALUES ($(join(map(x->join(topos(x),","), keys(chunk_labels)), "),(")));")
	end
	
	MySQL.execute!(mysql_conn, "START TRANSACTION;")
	for (cid, mappings) in chunk_labels
		if length(mappings) > 0
			x, y, z = topos(cid)
			res = MySQL.query(mysql_conn, "SELECT id FROM chunks WHERE x = $x AND y = $y AND z = $z;")
			MySQL.execute!(mysql_conn, "INSERT IGNORE INTO mappings (chunk_id, regiongraph_sv, chunkedgraph_sv) VALUES ($(res[:id][1]),$(join([join(pair, ',') for pair in mappings], "),($(res[:id][1]),")))")
		end
	end
	MySQL.execute!(mysql_conn, "COMMIT;")
end

function chunked_labelling(chunk::RelabeledChunk{T}) where T <: Integer
	# load existing labels of neighboring + center chunk
	print("Load labels from DB: "); tic(); 
	chunk.rg_to_cg_boundary = load_labels(chunk.watershed_original)
	chunk.rg_to_cg_complete = copy(chunk.rg_to_cg_boundary)
	cg_to_rg = Set{Label}()
	for c in values(chunk.rg_to_cg_complete)
		union!(cg_to_rg, values(c))
	end
	println(toq())

	print("Relabeling chunk: "); tic();
	chunk.watershed_relabeled = OffsetArray(Label, indices(chunk.watershed_original)...)
	nextsegid = Dict{ChunkID, UInt32}();
	for k in indices(chunk.watershed_original, 3), j in indices(chunk.watershed_original, 2), i in indices(chunk.watershed_original, 1)
		chunkid = ChunkedGraphs.world_to_chunk(i, j, k)
		if !haskey(chunk.rg_to_cg_complete, chunkid) # TODO: Create this outside
			chunk.rg_to_cg_complete[chunkid] = Dict{T, Label}()
			chunk.rg_to_cg_boundary[chunkid] = Dict{T, Label}()
		end
		if !haskey(nextsegid, chunkid)
			nextsegid[chunkid] = 1
		end

		# Don't relabel cell boundary
		if chunk.watershed_original[i, j, k] == 0
			chunk.watershed_relabeled[i, j, k] = 0
			continue
		end

		# Regiongraph to Chunked Graph
		if haskey(chunk.rg_to_cg_complete[chunkid], chunk.watershed_original[i, j, k])
			chunk.watershed_relabeled[i, j, k] = chunk.rg_to_cg_complete[chunkid][chunk.watershed_original[i, j, k]]
		else
			while tolabel(chunkid, UInt32(nextsegid[chunkid])) in cg_to_rg
				nextsegid[chunkid] += 1
			end

			chunk.watershed_relabeled[i, j, k] = tolabel(chunkid, UInt32(nextsegid[chunkid]))
			nextsegid[chunkid] += 1
			push!(cg_to_rg, chunk.watershed_relabeled[i, j, k])
			chunk.rg_to_cg_complete[chunkid][chunk.watershed_original[i, j, k]] = chunk.watershed_relabeled[i, j, k]
			if any((i, j, k) .% CG_CHUNKSIZE .== 0)
				chunk.rg_to_cg_boundary[chunkid][chunk.watershed_original[i, j, k]] = chunk.watershed_relabeled[i, j, k]
			end
		end
	end
	println(toq())

	# save new labels of neighboring + center chunk
	print("Storing labels to DB: "); tic(); 
	save_labels(chunk.rg_to_cg_boundary);
	println(toq())

	return chunk
end

function compute_regiongraph(watershed::AbstractArray{S, 3}, relabeled::AbstractArray{Label, 3}, agglomeration::AbstractArray{T, 3},
		regiongraph_edges::Set{AtomicEdge}) where {S <: Integer, T <: Integer}
	edges = Set{AtomicEdge}()
	@inbounds @simd for k in drop_last(indices(relabeled, 3))
		for j in drop_last(indices(relabeled, 2)), i in drop_last(indices(relabeled, 1))
			if relabeled[i, j, k] !== relabeled[i, j, k + 1] && agglomeration[i, j, k] === agglomeration[i, j, k + 1]
				if watershed[i, j, k] === watershed[i, j, k + 1]
					push!(edges, AtomicEdge(relabeled[i, j, k], relabeled[i, j, k + 1], INF_CAPACITY))
				else
					edge = getkey(regiongraph_edges.dict, AtomicEdge(watershed[i, j, k], watershed[i, j, k + 1]), nothing)
					if edge isa AtomicEdge
						push!(edges, AtomicEdge(relabeled[i, j, k], relabeled[i, j, k + 1], edge.affinity))
					end
				end
			end
		end
	end
	@inbounds @simd for k in drop_last(indices(relabeled, 3))
		for j in drop_last(indices(relabeled, 2)), i in drop_last(indices(relabeled, 1))
			if relabeled[i, j, k] !== relabeled[i, j + 1, k] && agglomeration[i, j, k] === agglomeration[i, j + 1, k]
				if watershed[i, j, k] === watershed[i, j + 1, k]
					push!(edges, AtomicEdge(relabeled[i, j, k], relabeled[i, j + 1, k], INF_CAPACITY))
				else
					edge = getkey(regiongraph_edges.dict, AtomicEdge(watershed[i, j, k], watershed[i, j + 1, k]), nothing)
					if edge isa AtomicEdge
						push!(edges, AtomicEdge(relabeled[i, j, k], relabeled[i, j + 1, k], edge.affinity))
					end
				end
			end
		end
	end
	@inbounds @simd for k in drop_last(indices(relabeled, 3))
		for j in drop_last(indices(relabeled, 2)), i in drop_last(indices(relabeled, 1))
			if relabeled[i, j, k] !== relabeled[i + 1, j, k] && agglomeration[i, j, k] === agglomeration[i + 1, j, k]
				if watershed[i, j, k] === watershed[i + 1, j, k]
					push!(edges, AtomicEdge(relabeled[i, j, k], relabeled[i + 1, j, k], INF_CAPACITY))
				else
					edge = getkey(regiongraph_edges.dict, AtomicEdge(watershed[i, j, k], watershed[i + 1, j, k]), nothing)
					if edge isa AtomicEdge
						push!(edges, AtomicEdge(relabeled[i, j, k], relabeled[i + 1, j, k], edge.affinity))
					end
				end
			end
		end
	end
	return collect(edges)
end

# struct RanStruct
# 	segA1::UInt64
# 	segB1::UInt64
# 	sum_aff1::Float32
# 	sum_area1::UInt64
# 	segA2::UInt64
# 	segB2::UInt64
# 	sum_aff2::Float32
# 	sum_area2::UInt64
# end
#
# The first 4 values are currently the same as the last 4 values...
# Affinity between segA1 and segB1 is sum_aff1 / sum_area1
# Values were written without padding, which is why we can't read them in Julia as Array

@memoize Dict function get_chunk_affinities(storage::StorageWrapper, chunk_path::AbstractString)
	edges = Vector{AtomicEdge}()

	f = storage.val[:get_file](chunk_path)
	if f isa Void
		warn("$chunk_path doesn't exist")
		return edges
	end

	buf = ZstdDecompressorStream(IOBuffer(f))
	sizehint!(edges, 500000)
	while !eof(buf)
		segA, segB = read(buf, UInt64, 2); aff = read(buf, Float32); area = read(buf, UInt64); skip(buf, 28);
		push!(edges, AtomicEdge(segA, segB, aff / area))
	end

	return edges
end

function get_chunkhierarchy_affinities(regiongraph::RegionGraphWrapper, slices::Tuple{UnitRange{T}, UnitRange{T}, UnitRange{T}}) where T
	chunk_range = map(x -> div.(x .- regiongraph.offset, regiongraph.chunksize), (minimum.(slices), maximum.(slices)))

	edges = Vector{AtomicEdge}()
	sizehint!(edges, 1000000)
	for l in 0:regiongraph.maxlevel
		for x in chunk_range[1][1]:chunk_range[2][1], y in chunk_range[1][2]:chunk_range[2][2], z in chunk_range[1][3]:chunk_range[2][3]
			chunk_path = "complete_edges_$(l)_$(x)_$(y)_$(z).data.zst"
			append!(edges, get_chunk_affinities(regiongraph.storage, chunk_path))
		end
		chunk_range = (div.(first(chunk_range), 2), div.(last(chunk_range), 2))
	end

	return edges
end

function edge_task(watershed::CloudVolumeWrapper, agglomeration::CloudVolumeWrapper, regiongraph::RegionGraphWrapper,
		output_storage::StorageWrapper, slices::Tuple{UnitRange{T}, UnitRange{T}, UnitRange{T}}) where T<:Integer

	print("Get regiongraph edges: "); tic();
	regiongraph_edges = get_chunkhierarchy_affinities(regiongraph, map(drop_last, slices))
	println(toq())

	print("Converting regiongraph edges to Set: "); tic();
	regiongraph_edges = Set(regiongraph_edges)
	println(toq())

	print("Get cutouts: "); tic();
	watershed_cutout = OffsetArray(watershed[slices...], slices...)
	agglomeration_cutout = OffsetArray(agglomeration[slices...], slices...)
	println(toq())

	chunk = RelabeledChunk(watershed_cutout)
	chunk = chunked_labelling(chunk)

	print("Compute Regiongraph: "); tic();
	edges = compute_regiongraph(watershed_cutout, chunk.watershed_relabeled, agglomeration_cutout, regiongraph_edges)
	println(toq())

	print("Write results: "); tic();
	chunkid = ChunkedGraphs.world_to_chunk(map(first, slices)...)
	edge_buf = IOBuffer()
	rg_to_cg_buf = IOBuffer()
	write(edge_buf, edges)
	write(rg_to_cg_buf, sort(collect(chunk.rg_to_cg_complete[chunkid])))

	output_storage.val[:put_file](file_path="$(stringify(slices))_atomicedges.bin", content = PyCall.pybytes(edge_buf.data[1:edge_buf.size]))
	output_storage.val[:put_file](file_path="$(stringify(slices))_rg2cg.bin", content = PyCall.pybytes(rg_to_cg_buf.data[1:rg_to_cg_buf.size]))
	output_storage.val[:wait]()
	println(toq())
end

function edge_task_bundle(watershed_path::AbstractString, agglomeration_path::AbstractString,
		regiongraph_path::AbstractString, output_path::AbstractString,
		slices::Tuple{UnitRange{T}, UnitRange{T}, UnitRange{T}}) where T <: Integer
	slices = map(drop_last, slices) # julia indexing is inclusive

	watershed = CloudVolumeWrapper(watershed_path; bounded=false)
	agglomeration = CloudVolumeWrapper(agglomeration_path; bounded=false)

	regiongraph = RegionGraphWrapper(
		StorageWrapper(regiongraph_path),
		RG_CHUNKSIZE,
		(offset(watershed)...),
		Int(ceil(log2(maximum( div.(size(watershed)[1:3], RG_CHUNKSIZE) )))) # max octree level
	)
	output_storage = StorageWrapper(output_path)

	for x_start in first(slices[1]) : CG_CHUNKSIZE[1] : last(slices[1]) - 1
		for y_start in first(slices[2]) : CG_CHUNKSIZE[2] : last(slices[2]) - 1
			for z_start in first(slices[3]) : CG_CHUNKSIZE[3] : last(slices[3]) - 1
				subslice = (x_start : x_start + CG_CHUNKSIZE[1],
				            y_start : y_start + CG_CHUNKSIZE[2],
				            z_start : z_start + CG_CHUNKSIZE[3])
				edge_task(watershed, agglomeration, regiongraph, output_storage, subslice)
			end
		end
	end
end

watershed_path, agglomeration_path, regiongraph_path, output_path, slice_str = ARGS

# watershed_path = "gs://neuroglancer/ranl/pinky40_watershed_test"
# agglomeration_path = "gs://neuroglancer/ranl/pinky40_agglomeration_test"
# regiongraph_path = "gs://ranl/pinky40_agglomeration_test/region_graph"
# output_path = "gs://nkem/pinky40_agglomeration_test/region_graph"

# slice_str = "10240-11265_7680-8705_0-1025"

edge_task_bundle(watershed_path, agglomeration_path, regiongraph_path, output_path, toslice(slice_str))
