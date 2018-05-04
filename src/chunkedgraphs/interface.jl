function children(cgraph::ChunkedGraph, vertices::Vector{Vertex})
	reqlength = mapreduce(v -> length(v.children), +, 0, vertices)
	if reqlength == 0
		return Label[]
	end
	
	vertexlabels = Vector{Label}(reqlength)
	i = 1
	for v in vertices
		vertexlabels[i : i + length(v.children) - 1] = v.children
		i += length(v.children)
	end
	return vertexlabels
end

function leaves!(cgraph::ChunkedGraph, vertex::Vertex, stop_lvl::Integer = 1, bbox::Union{Cuboid, Void} = nothing)
	@assert tolevel(tochunkid(vertex)) >= stop_lvl

	if tolevel(tochunkid(vertex)) == stop_lvl
		return [vertex.label]
	end

	vertices = [vertex]
	lvl = tolevel(tochunkid(vertex))
	while lvl > stop_lvl + 1
		vertexlabels = children(cgraph, vertices)
		if bbox !== nothing
			filter!(v -> overlaps(tocuboid(cgraph, tochunkid(v)), bbox::Cuboid), vertexlabels)
		end
		vertices = map(lbl -> getvertex!(cgraph, lbl), vertexlabels)
		lvl -= 1
	end
	return children(cgraph, vertices)
end

function promote!(cgraph::ChunkedGraph, vertex::Vertex)
	c = getchunk!(cgraph, tochunkid(vertex))
	@assert tolevel(c) < cgraph.MAX_DEPTH

	@assert c.clean
	@assert vertex.parent == NULL_LABEL
	@assert length(incident_edges(c.graph, vertex.label)) == 0

	l = uniquelabel!(c.parent)
	pv = Vertex(l, NULL_LABEL, Label[vertex.label])
	vertex.parent = pv.label

	@assert tochunkid(pv) == c.parent.id
	add_vertex!(c.parent.graph, pv.label)
	c.parent.vertices[pv.label] = pv
	c.parent.modified = true

	return pv
end

function promote_to_lca!(cgraph::ChunkedGraph, vertex1::Vertex, vertex2::Vertex)
	if tochunkid(vertex1) == tochunkid(vertex2)
		return (vertex1, vertex2)
	else
		if vertex1.parent == NULL_LABEL
			promote!(cgraph, vertex1)
		end
		if vertex2.parent == NULL_LABEL
			promote!(cgraph, vertex2)
		end
		parent1 = getvertex!(cgraph, vertex1.parent)
		parent2 = getvertex!(cgraph, vertex2.parent)
		return promote_to_lca!(cgraph, parent1, parent2)
	end
end

function root!(cgraph::ChunkedGraph, vertex::Vertex)
	if isroot(vertex)
		return vertex
	else
		return root!(cgraph, getvertex!(cgraph, vertex.parent))
	end
end

n_processed = 0
function update!(cgraph::ChunkedGraph)
	global n_processed
	n_processed = 0
	gc_enable(false)
	update!(getchunk!(cgraph, cgraph.TOP_ID))
	gc_enable(true)
end

function add_atomic_vertex!(cgraph::ChunkedGraph, lbl::Label)
	@assert tolevel(tochunkid(lbl)) == 1 "Vertex label at level $(tolevel(tochunkid(label))), expected 1."

	c = getchunk!(cgraph, tochunkid(lbl))
	if haskey(c.vertices, lbl)
		#TODO: warn user
		return
	end

	v = Vertex(lbl, NULL_LABEL, EMPTY_LABEL_LIST)
	push!(c.added_vertices, v)
	touch!(c)
end

function add_atomic_vertices!(cgraph::ChunkedGraph, lbls::Vector{Label})
	gc_enable(false)
	for lbl in lbls
		add_atomic_vertex!(cgraph, lbl)
	end
	gc_enable(true)
end

function add_atomic_edge!(cgraph::ChunkedGraph, edge::AtomicEdge)
	c = getchunk!(cgraph, lca(cgraph, tochunkid(edge.u), tochunkid(edge.v)))
	push!(c.added_edges, edge)
	touch!(c)
end

function add_atomic_edges!(cgraph::ChunkedGraph, edges::Vector{AtomicEdge})
	gc_enable(false)
	for edge in edges
		add_atomic_edge!(cgraph, edge)
	end
	gc_enable(true)
end

function delete_atomic_edge!(cgraph::ChunkedGraph, edge::AtomicEdge)
	c = getchunk!(cgraph, lca(cgraph, tochunkid(edge.u), tochunkid(edge.v)))
	push!(c.deleted_edges, edge)
	touch!(c)
end
