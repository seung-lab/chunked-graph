# TODO: Proper structure, more tests, add Coverage
push!(LOAD_PATH, dirname(@__FILE__))
include("../src/chunkedgraphs/ChunkedGraphs.jl")
using ChunkedGraphs
using Base.Test

function test_cases()
	@testset "all_tests" begin
		if !isdir("/tmp/graph")
			mkdir("/tmp/graph")
		end

		@testset "add test_add_atomic_node" begin
			G = ChunkedGraph("/tmp/graph")
			label = tolabel(1,0,0,0,1)
			add_atomic_vertex!(G, label)
			update!(G)
			@test getvertex!(G, label) == Vertex(label, NULL_LABEL, EMPTY_LABEL_LIST)

			@test_throws KeyError getvertex!(G, tolabel(1,0,0,0,2))
		end

		@testset "test_circle" begin
			G = ChunkedGraph("/tmp/graph")
			add_atomic_vertex!(G, tolabel(1,0,0,0,3))
			add_atomic_vertex!(G, tolabel(1,0,0,0,2))
			add_atomic_vertex!(G, tolabel(1,0,0,0,1))
			update!(G)

			add_atomic_edge!(G, AtomicEdge(tolabel(1,0,0,0,1), tolabel(1,0,0,0,2), 1.f0))
			add_atomic_edge!(G, AtomicEdge(tolabel(1,0,0,0,1), tolabel(1,0,0,0,3), 1.f0))
			add_atomic_edge!(G, AtomicEdge(tolabel(1,0,0,0,2), tolabel(1,0,0,0,3), 1.f0))
			update!(G)

			@test getvertex!(G, tolabel(1,0,0,0,1)).parent == getvertex!(G, tolabel(1,0,0,0,2)).parent
			@test getvertex!(G, tolabel(1,0,0,0,2)).parent == getvertex!(G, tolabel(1,0,0,0,3)).parent
		end

		@testset "test_circle_external_edge" begin
			G = ChunkedGraph("/tmp/graph")
			add_atomic_vertex!(G, tolabel(1,0,0,1,3))
			add_atomic_vertex!(G, tolabel(1,0,0,0,2))
			add_atomic_vertex!(G, tolabel(1,0,0,0,1))

			add_atomic_edge!(G, AtomicEdge(tolabel(1,0,0,1,3), tolabel(1,0,0,0,2), 1.f0))
			update!(G)
			add_atomic_edge!(G, AtomicEdge(tolabel(1,0,0,0,2), tolabel(1,0,0,0,1), 1.f0))
			update!(G)
			@test root!(G, getvertex!(G, tolabel(1,0,0,1,3))) == root!(G, getvertex!(G, tolabel(1,0,0,0,2)))
			@test root!(G, getvertex!(G, tolabel(1,0,0,1,3))) == root!(G, getvertex!(G, tolabel(1,0,0,0,1)))
		end

		@testset "delete_edge_same_chunk" begin
			G = ChunkedGraph("/tmp/graph")
			add_atomic_vertex!(G, tolabel(1,0,0,0, 1) )
			add_atomic_vertex!(G, tolabel(1,0,0,0, 2) )
			add_atomic_edge!(G, AtomicEdge(tolabel(1,0,0,0,1), tolabel(1,0,0,0,2), 1.f0))
			update!(G)
			@test root!(G, getvertex!(G, tolabel(1,0,0,0,1))) == root!(G, getvertex!(G, tolabel(1,0,0,0,2)))

			delete_atomic_edge!(G, AtomicEdge(tolabel(1,0,0,0,1), tolabel(1,0,0,0,2)))
			update!(G)
			@test root!(G, getvertex!(G, tolabel(1,0,0,0,1))) != root!(G, getvertex!(G, tolabel(1,0,0,0,2)))
		end

		@testset "delete_edge_different_chunk" begin
			G = ChunkedGraph("/tmp/graph")
			u = tolabel(1,0,0,0,1)
			v = tolabel(1,0,0,1,2)
			add_atomic_vertex!(G, u)
			add_atomic_vertex!(G, v)
			add_atomic_edge!(G, AtomicEdge(u,v, 1.f0))
			update!(G)
			@test root!(G, getvertex!(G,u)) == root!(G, getvertex!(G,v))

			delete_atomic_edge!(G, AtomicEdge(u,v))
			update!(G)
			@test root!(G, getvertex!(G,u)) != root!(G, getvertex!(G,v))
			@test length(root!(G, getvertex!(G, u)).children) == 1
			@test length(root!(G, getvertex!(G, v)).children) == 1
		end

		@testset "test_3_node_delete" begin
			G = ChunkedGraph("/tmp/graph")
			add_atomic_vertex!(G, tolabel(1,0,0,0,1) )
			add_atomic_vertex!(G, tolabel(1,0,0,1,2) )
			add_atomic_vertex!(G, tolabel(1,0,0,3,3) )

			add_atomic_edge!(G, AtomicEdge(tolabel(1,0,0,0,1), tolabel(1,0,0,1,2), 1.f0))
			add_atomic_edge!(G, AtomicEdge(tolabel(1,0,0,1,2), tolabel(1,0,0,3,3), 1.f0))
			update!(G)

			delete_atomic_edge!(G, AtomicEdge(tolabel(1,0,0,0,1), tolabel(1,0,0,1,2)))
			update!(G)
			@test root!(G, getvertex!(G, tolabel(1,0,0,0,1))) != root!(G, getvertex!(G, tolabel(1,0,0,1,2)))
			@test root!(G, getvertex!(G, tolabel(1,0,0,1,2))) == root!(G, getvertex!(G, tolabel(1,0,0,3,3)))

			@test length(root!(G, getvertex!(G, tolabel(1,0,0,0,1))).children) == 1
		end


		@testset "two_node_mincut!" begin
			G = ChunkedGraph("/tmp/graph")
			u = tolabel(1,0,0,0,1)
			v = tolabel(1,0,0,0,2)

			add_atomic_vertex!(G, u )
			add_atomic_vertex!(G, v )

			add_atomic_edge!(G, AtomicEdge(u, v, 1.f0))

			update!(G)
			@test Set(mincut!(G, u, v)) == Set([AtomicEdge(u,v)])

			println(mincut!(G, u, u))
			# @test_throws KeyError mincut!(G, v, v)
			@test_throws KeyError mincut!(G, u, tolabel(1,0,0,0,3))
			@test_throws KeyError mincut!(G, u, tolabel(1,0,0,5,1))
		end


		@testset "triangle_mincut!" begin
			G = ChunkedGraph("/tmp/graph")
			u = tolabel(1,0,0,0,1)
			v = tolabel(1,0,0,0,2)
			w = tolabel(1,0,0,0,3)

			add_atomic_vertex!(G, u )
			add_atomic_vertex!(G, v )
			add_atomic_vertex!(G, w )

			add_atomic_edge!(G, AtomicEdge(u, v, 1.f0))
			add_atomic_edge!(G, AtomicEdge(u, w, 1.f0))
			add_atomic_edge!(G, AtomicEdge(v, w, 1.f0))

			update!(G)

			@test Set(mincut!(G, u, v)) == Set([AtomicEdge(u,v),AtomicEdge(v,w)])
			@test Set(mincut!(G, [u,w], [v])) == Set([AtomicEdge(u,v),AtomicEdge(v,w)])
			@test Set(mincut!(G, [u], [v,w])) == Set([AtomicEdge(u,v),AtomicEdge(u,w)])
		end

		@testset "chunk_mincut!" begin
			G = ChunkedGraph("/tmp/graph")
			u = tolabel(1,0,0,0,1)
			v = tolabel(1,0,0,0,2)
			w = tolabel(1,0,0,1,3)

			add_atomic_vertex!(G, u )
			add_atomic_vertex!(G, v )
			add_atomic_vertex!(G, w )

			add_atomic_edge!(G, AtomicEdge(u, v, 1.f0))
			add_atomic_edge!(G, AtomicEdge(u, w, 1.f0))
			add_atomic_edge!(G, AtomicEdge(v, w, 1.f0))

			update!(G)

			@test Set(mincut!(G, u, v)) == Set([AtomicEdge(u,v),AtomicEdge(v,w)])
			@test Set(mincut!(G, [u,w], [v])) == Set([AtomicEdge(u,v),AtomicEdge(v,w)])
			@test Set(mincut!(G, [u], [v,w])) == Set([AtomicEdge(u,v),AtomicEdge(u,w)])
		end

		@testset "affinity_mincut!" begin
			G = ChunkedGraph("/tmp/graph")
			source = tolabel(1,0,0,0,1)
			a1 = tolabel(1,0,0,0,2)
			a2 = tolabel(1,0,0,0,3)
			a3 = tolabel(1,0,0,0,4)
			b = tolabel(1,0,0,0,5)
			sink = tolabel(1,0,0,0,6)

			add_atomic_vertex!(G, source)
			add_atomic_vertex!(G, sink)
			add_atomic_vertex!(G, a1)
			add_atomic_vertex!(G, a2)
			add_atomic_vertex!(G, a3)
			add_atomic_vertex!(G, b)

			add_atomic_edge!(G, AtomicEdge(source, a1, 0.25f0))
			add_atomic_edge!(G, AtomicEdge(source, a2, 0.25f0))
			add_atomic_edge!(G, AtomicEdge(source, a3, 0.25f0))

			add_atomic_edge!(G, AtomicEdge(a1, a2, 1.0f0))
			add_atomic_edge!(G, AtomicEdge(a2, a3, 1.0f0))
			add_atomic_edge!(G, AtomicEdge(a3, a1, 1.0f0))

			add_atomic_edge!(G, AtomicEdge(a1, b, 1.0f0))
			add_atomic_edge!(G, AtomicEdge(a2, b, 1.0f0))
			add_atomic_edge!(G, AtomicEdge(a3, b, 1.0f0))

			add_atomic_edge!(G, AtomicEdge(b, sink, 0.76))
			update!(G)

			@test Set(mincut!(G, source, sink)) == Set([AtomicEdge(source,a1),AtomicEdge(source,a2),AtomicEdge(source,a3)])
		end

		@testset "multi_split" begin
			# Two triangles connected over a small bridge (x1-x2)
			G = ChunkedGraph("/tmp/graph")
			a1 = tolabel(1,0,0,0,1)
			a2 = tolabel(1,0,0,0,2)
			x1 = tolabel(1,0,0,0,3)
			x2 = tolabel(1,0,0,0,4)
			b1 = tolabel(1,0,0,0,5)
			b2 = tolabel(1,0,0,0,6)

			add_atomic_vertex!(G, a1 )
			add_atomic_vertex!(G, a2 )
			add_atomic_vertex!(G, x1 )
			add_atomic_vertex!(G, x2 )
			add_atomic_vertex!(G, b1 )
			add_atomic_vertex!(G, b2 )

			add_atomic_edge!(G, AtomicEdge(a1, a2, 1.f0))
			add_atomic_edge!(G, AtomicEdge(a1, x1, 1.f0))
			add_atomic_edge!(G, AtomicEdge(a2, x1, 1.f0))
			add_atomic_edge!(G, AtomicEdge(x1, x2, 1.f0))
			add_atomic_edge!(G, AtomicEdge(x2, b1, 1.f0))
			add_atomic_edge!(G, AtomicEdge(x2, b2, 1.f0))
			add_atomic_edge!(G, AtomicEdge(b1, b2, 1.f0))

			update!(G)

			@test Set(mincut!(G, [a1,a2], [b1])) == Set([AtomicEdge(x1,x2)])
			@test Set(mincut!(G, [a1,b1], [b2])) == Set([AtomicEdge(x2,b2),AtomicEdge(b1,b2)])
			@test Set(mincut!(G, [a1,b1], [a2,b2])) == Set([AtomicEdge(a1,a2),AtomicEdge(b1,b2),AtomicEdge(a2,x1),AtomicEdge(x2,b2)])
		end

		@testset "multi_split_omni" begin
			#= Copied from Omni
			* Reference trinary tree diagram for the following tests
			* Labelled Using vertex INDEX NOT segID
			*                                         1
			*                                         ^
			*               /                         |                         \
			*              2                          3                          4
			*              ^                          ^                          ^
			*      /       |       \          /       |       \          /       |       \
			*     5        6        7        8        9       10        11      12       13
			*     ^        ^        ^        ^        ^        ^        ^        ^        ^
			*   / |  \   / |  \   / |  \   / |  \   / |  \   / |  \   / |  \   / |  \   / |  \
			*  14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
			=#
			G = ChunkedGraph("/tmp/graph")
			v = map(x->tolabel(1,0,0,0,x), 1:40)
			for vertex in v
				add_atomic_vertex!(G, vertex)
			end

			for i=1:13
				add_atomic_edge!(G, AtomicEdge(v[i], v[3*i-1], 0.5f0))
				add_atomic_edge!(G, AtomicEdge(v[i], v[3*i], 0.5f0))
				add_atomic_edge!(G, AtomicEdge(v[i], v[3*i+1], 0.5f0))
			end

			update!(G)

			@test Set(mincut!(G, [v[17], v[6], v[3]], [v[18], v[4], v[12]])) == Set([AtomicEdge(v[1],v[4]),AtomicEdge(v[6], v[18])])
		end

		@testset "supervoxels_not_splitted" begin
			#=
			Currently are datasets have unique
			supervoxels
			this means seg_id (the first 32bits of a label)
			are unique.

			We will need to remove this limitation eventually
			but for now, we will make sure we won't split 
			supervoxels with same seg id
			=#

			G = ChunkedGraph("/tmp/graph")
			u = tolabel(1,0,0,0,1)
			v = tolabel(1,0,0,1,1)
			add_atomic_vertex!(G, u )
			add_atomic_vertex!(G, v )
			add_atomic_edge!(G, AtomicEdge(u, v, 1.f0))
			update!(G)

			@test isempty(mincut!(G, u, v))
		end
	end

end

test_cases()
