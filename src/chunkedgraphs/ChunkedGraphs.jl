__precompile__(true)
push!(LOAD_PATH, dirname(@__FILE__))
module ChunkedGraphs

export ChunkedGraph, Chunk
export AtomicEdge, Vertex
export Label, Affinity

export INF_CAPACITY, NULL_LABEL, EMPTY_LABEL_LIST

export tolabel, topos, tosegid, tolevel, tochunkid
export add_atomic_vertex!, add_atomic_vertices!, add_atomic_edge!, add_atomic_edges!, delete_atomic_edge!
export mincut!, update!, getvertex!, getchunk!, root!, leaves!, save!
export getsupervoxelat

include("core.jl")
include("utils.jl")
include("vertex.jl")
include("edge.jl")
include("multigraph.jl")
include("chunk.jl")
include("io.jl")
include("interface.jl")
include("mincut.jl")
include("segmentpicker.jl")
end
