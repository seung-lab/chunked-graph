using DataStructures

const Label            = UInt64
const Affinity         = Float32

const ChunkID          = UInt32
const SegmentID        = UInt32

const Cuboid           = Tuple{UnitRange{Int}, UnitRange{Int}, UnitRange{Int}}

const INF_CAPACITY     = typemax(Affinity)
const NULL_LABEL       = typemax(Label)
const EMPTY_LABEL_LIST = Vector{Label}()

# TODO: level, x, y, z are currently fixed to 8 bit, respectively. Should be adjustable.
const low_mask_8       = UInt32(0x000000FF)
const low_mask_32      = UInt64(0x00000000FFFFFFFF)

@inline function tolabel(chk::ChunkID, seg::SegmentID)
	return Label(seg) | (Label(chk) << 32)
end

@inline function tolabel(lvl::Integer, x::Integer, y::Integer, z::Integer, seg::Integer)
	return tolabel(tochunkid(lvl, x, y, z), SegmentID(seg))
end

@inline function tochunkid(lvl::UInt32, x::UInt32, y::UInt32, z::UInt32)
	return ChunkID((lvl << 24) | (x << 16) | (y << 8) | z)
end

@inline function tochunkid(lvl::Integer, x::Integer, y::Integer, z::Integer)
	return tochunkid(UInt32(lvl), UInt32(x), UInt32(y), UInt32(z))
end

@inline function tochunkid(lbl::Label)
	return ChunkID(lbl >> 32)
end

@inline function tosegid(lbl::Label)
	return SegmentID(lbl & low_mask_32)
end

@inline function tolevel(chk::ChunkID)
	return UInt8(chk >> 24)
end

@inline function tolevel(lbl::Label)
	return UInt8(lbl >> 56)
end

@inline function topos(chk::ChunkID)
	return UInt8((chk >> 16) & low_mask_8), UInt8((chk >> 8) & low_mask_8), UInt8(chk & low_mask_8)
end

@inline function topos(lbl::Label)
	return UInt8((lbl >> 48) & low_mask_8), UInt8((lbl >> 40) & low_mask_8), UInt8((lbl >> 32) & low_mask_8)
end

"Creates a bounding box as `Tuple{UnitRange{Int}, UnitRange{Int}, UnitRange{Int}}`. Coordinates are *chunk* coordinates."
function tocuboid(chk::ChunkID)
	@assert tolevel(chk) >= 1
	if chk === TOP_ID || chk === SECOND_ID
		return (typemin(Int):typemax(Int), typemin(Int):typemax(Int), typemin(Int):typemax(Int))::Cuboid
	else
		mult = 2^(tolevel(chk) - 1)
		x, y, z = topos(chk)
		return (x * mult : (x + 1) * mult, y * mult : (y + 1) * mult, z * mult : (z + 1) * mult)::Cuboid
	end
end

# TODO: level > 1, switch to BoundingBoxes.jl?
"Creates a bounding box as `Tuple{UnitRange{Int}, UnitRange{Int}, UnitRange{Int}}`. Coordinates are *chunk* coordinates."
function tocuboid(lbls::Vector{Label}, dilate::Int = 0)
	min_x, min_y, min_z = typemax(Int), typemax(Int), typemax(Int)
	max_x, max_y, max_z = 0, 0, 0

	for lbl in lbls
		@assert tolevel(tochunkid(lbl)) == 1
		x, y, z = topos(tochunkid(lbl))
		min_x = min(min_x, x); max_x = max(max_x, x)
		min_y = min(min_y, y); max_y = max(max_y, y)
		min_z = min(min_z, z); max_z = max(max_z, z)
	end

	return (min_x - dilate : max_x + dilate,
			min_y - dilate : max_y + dilate,
			min_z - dilate : max_z + dilate)::Cuboid
end

@inline function overlaps(r1::UnitRange, r2::UnitRange)
	return r1.start <= r2.stop && r2.start <= r1.stop
end

@inline function overlaps(c1::Cuboid, c2::Cuboid)
	@inbounds return overlaps(c1[1], c2[1]) && overlaps(c1[2], c2[2]) && overlaps(c1[3], c2[3])
end

@inline function isroot(chunkid::ChunkID)
	return chunkid === TOP_ID
end

"Calculates the parent's ChunkID for a given chunk ID"
function parent(chunkid::ChunkID)
	if tolevel(chunkid) >= MAX_DEPTH
		return TOP_ID
	elseif tolevel(chunkid) == MAX_DEPTH - 1
		return SECOND_ID
	else
		x, y, z = topos(chunkid)
		return tochunkid(tolevel(chunkid) + 1, fld(x, 2), fld(y, 2), fld(z, 2))
	end
end

"Calculates the last/lowes common ancestor's ChunkID for two given chunk IDs"
function lca(chunkid1::ChunkID, chunkid2::ChunkID)
	# TODO: Ensure chunks are on same level
	@assert tolevel(chunkid1) === tolevel(chunkid2)
	if chunkid1 === chunkid2
		return chunkid1
	else
		return lca(parent(chunkid1), parent(chunkid2))
	end
end

"Create a string of form 'lvl_x_y_z' for a given ChunkID"
function stringify(chunkid::ChunkID)
	x, y, z = topos(chunkid)
	return String("$(tolevel(chunkid))_$(x)_$(y)_$(z)")
end

@inline function world_to_chunk(x::Integer, y::Integer, z::Integer)
	return tochunkid(1, fld(x, CHUNK_SIZE[1]), fld(y, CHUNK_SIZE[2]), fld(z, CHUNK_SIZE[3]))
end


const MAX_DEPTH        = 8
const TOP_ID           = tochunkid(MAX_DEPTH + 1, 0, 0, 0)
const SECOND_ID        = tochunkid(MAX_DEPTH, 0, 0, 0)

const CHUNK_SIZE       = (512, 512, 64)

const CACHESIZE = 40000
eviction_mode = false #TODO: fix this

mutable struct ChunkedGraph{C} # {C} is necessary until Julia supports forward declaration of Chunk
	chunks::Dict{ChunkID, C}
	lastused::PriorityQueue{ChunkID, Float64}
	path::AbstractString
end

function ChunkedGraph(graphpath::AbstractString)
	@assert isdir(graphpath)
	return ChunkedGraph{Chunk}(
		Dict{ChunkID, Chunk}(),
		PriorityQueue{ChunkID, Float64}(),
		graphpath
	)
end
