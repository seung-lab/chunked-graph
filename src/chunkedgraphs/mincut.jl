using LightGraphsFlows

function common_parent_vertices!(cgraph::ChunkedGraph, lbls::Vector{Label})
	@assert length(lbls) >= 1

	# Get unique set of root_vertices (not necessarily on top-level)
	root_vertex_set = Set{Vertex}()
	for lbl in lbls
		@assert tolevel(tochunkid(lbl)) == 1
		push!(root_vertex_set, root!(cgraph, getvertex!(cgraph, lbl)))
	end

	root_vertices = collect(root_vertex_set)

	# Make sure all root_vertices are on the same level
	max_level = max(map(r->tolevel(tochunkid(r.label)), root_vertices)...)
	for i in eachindex(root_vertices)
		while tolevel(tochunkid(root_vertices[i].label)) < max_level
			root_vertices[i] = promote!(cgraph, root_vertices[i])
		end
	end

	# Make sure all root_vertices are in the same chunk (already on the same level)
	while length(unique(map(r->tochunkid(r), root_vertices))) !== 1
		root_vertices = map(r->promote!(cgraph, r), root_vertices)
	end

	return root_vertices
end


function induced_subgraph!(cgraph::ChunkedGraph, chunk::Chunk, vertices::Vector{Label}, bbox::Cuboid)
	atomicedges = AtomicEdge[]
	chunks = [chunk]
	lvl = tolevel(first(chunks))

	while lvl > 0
		# Add all induced edges of all important vertices in all chunks on the current level to the existing list of atomicedges
		#append!(atomicedges, vcat([vcat(values(induced_edges(c.graph, vertices))...) for c in chunks]...)); # arrays
		for c in chunks
			append!(atomicedges, induced_atomic_edges(c.graph, vertices))
		end

		if lvl > 2
			# From the current set of vertices, collect all child vertices
			vertices = convert(Vector{Label}, vcat([vcat([c.vertices[v].children for v in filter(v->haskey(c.vertices, v), vertices)]...) for c in chunks]...));

			# Make sure we got all the necessary chunks in memory
			foreach(v->getchunk!(cgraph, tochunkid(v)), vertices)

			# From the current set of chunks, collect all child chunks that still lie within the ROI
			chunks = filter(subc->overlaps(tocuboid(tochunkid(subc)), bbox), vcat([c.children for c in chunks]...))
		end
		lvl -= 1
	end

	atomicvertices = Set(vertices)
	return collect(atomicvertices), filter(e->e.u in atomicvertices && e.v in atomicvertices, atomicedges)
end

function mincut!(cgraph::ChunkedGraph, source::Label, sink::Label)
	return mincut!(cgraph, [source], [sink])
end

function mincut!(cgraph::ChunkedGraph, sources::Vector{Label}, sinks::Vector{Label})
	bbox = tocuboid(vcat(sources, sinks))
	root_vertices = common_parent_vertices!(cgraph, vcat(sources, sinks))
	c = getchunk!(cgraph, parent(tochunkid(root_vertices[1].label)))
	atomic_vertices, atomic_edges = induced_subgraph!(cgraph, c, map(r->r.label, root_vertices), bbox)

	# Relabel atomic vertex labels to get a dense matrix for LightGraphs
	encode = Dict{Label, Int}()
	for (i, lbl) in enumerate(atomic_vertices)
		encode[lbl] = i
	end

	N = length(atomic_vertices)
	fake_source = N + 1
	fake_sink = N + 2
	N += 2
	flow_graph = LightGraphs.DiGraph(N)
	capacities = zeros(Affinity, (N, N))
	for atomic_edge in atomic_edges
		u = atomic_edge.u
		v = atomic_edge.v
		affinity = atomic_edge.affinity

		LightGraphs.add_edge!(flow_graph, encode[u], encode[v])
		LightGraphs.add_edge!(flow_graph, encode[v], encode[u])

		# Don't split supervoxels at chunk boundaries
		# FIXME: this hack only works because segment are currently unique
		# across the whole dataset
		capacities[encode[u], encode[v]] = tosegid(u) == tosegid(v) ? INF_CAPACITY : affinity 
		capacities[encode[v], encode[u]] = tosegid(u) == tosegid(v) ? INF_CAPACITY : affinity
	end

	# create a fake source and add edges with infinite weight to all sources
	for source in sources
		LightGraphs.add_edge!(flow_graph, encode[source], fake_source)
		LightGraphs.add_edge!(flow_graph, fake_source, encode[source])

		capacities[encode[source], fake_source] = INF_CAPACITY
		capacities[fake_source, encode[source]] = INF_CAPACITY
	end

	# create a fake sink
	for sink in sinks
		LightGraphs.add_edge!(flow_graph, encode[sink], fake_sink)
		LightGraphs.add_edge!(flow_graph, fake_sink, encode[sink])

		capacities[encode[sink], fake_sink] = INF_CAPACITY
		capacities[fake_sink, encode[sink]] = INF_CAPACITY
	end

	f, _, labels = LightGraphsFlows.multiroute_flow(
		flow_graph, fake_source, fake_sink, capacities,
		flow_algorithm = LightGraphsFlows.BoykovKolmogorovAlgorithm(),
		routes = 1)

	# No split found, or no split allowed
	if f < 0.00001 || f == INF_CAPACITY
		return AtomicEdge[]
	end

	# labels contains either 1, 2 or 0
	# 1 means that the vertex in position i should be part of the source
	# 2 means that it should be part of the sink
	# 0 means that it could be either 1 or 2
	# We will transform all 0 to 1, so that this function always
	# return two parts
	labels = map(x->x == 0 ? 1 : x, labels)
	edges_to_cut = filter(e->labels[encode[e.u]] != labels[encode[e.v]], atomic_edges)
	
	# DISCUSSION:
	# We connected sinks and sources each by using a fake source/sink.
	# It is possible that real sources/sinks become disconnected after
	# cutting these edges - or already were disconnected in the beginning.
	# Should we create edges between all sinks and sources, respectively?
	return edges_to_cut
end
