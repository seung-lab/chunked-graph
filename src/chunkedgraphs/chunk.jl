using IterTools
using FileLocks

mutable struct Chunk
	cgraph::ChunkedGraph{Chunk}
	id::ChunkID
	graph::MultiGraph{Label,AtomicEdge,CompositeEdgeSet}
	vertices::Dict{Label, Vertex}
	parent::Union{Chunk, Void}
	children::Vector{Chunk}
	added_vertices::Dict{Label,Vertex}
	deleted_vertices::Dict{Label,Vertex}
	added_edges::Set{AtomicEdge}
	deleted_edges::Set{AtomicEdge}
	clean::Bool
	modified::Bool
	max_label::SegmentID
	flock::FileLock

	#todo: lazy loading for vertices

	function Chunk(cgraph::ChunkedGraph, chunkid::ChunkID, vertices::Dict{Label,Vertex}, graph::MultiGraph, max_label::Label)
		prefix = stringify(chunkid)
		fl = FileLock(expanduser(joinpath(cgraph.path, String("$(prefix).lock"))))
		lock(fl)

		c = new(cgraph, chunkid, graph, vertices, nothing, Chunk[], Dict{Label,Vertex}(), Dict{Label,Vertex}(), Set{AtomicEdge}(), Set{AtomicEdge}(), true, false, max_label, fl)

		if !isroot(chunkid)
			par = parent!(c)
			push!(par.children, c)
			c.parent = par
		end
		return c
	end
end

@inline function add_vertex!(c::Chunk, v::Vertex)
	@assert !(v.label in keys(c.deleted_vertices))
	@assert parent(tochunkid(v.label)) == c.id
	c.added_vertices[v.label] = v
	touch!(c)
end
@inline function delete_vertex!(c::Chunk, v::Vertex)
	@assert haskey(c.vertices,v)
	@assert parent(tochunkid(v.label)) == c.id
	c.deleted_vertices[v.label] = v
	touch!(c)
end
@inline function delete_vertex!(c::Chunk, l::Label)
	if haskey(c.vertices,l)
		v=c.vertices[l]::Vertex
	elseif haskey(c.added_vertices,l)
		v=c.added_vertices[l]::Vertex
	else
		@assert false
	end
	@assert parent(tochunkid(v.label)) == c.id
	c.deleted_vertices[l] = v
	touch!(c)
end
@inline function add_edge!(c::Chunk, e::AtomicEdge)
	@assert !(e in c.deleted_edges)
	@assert begin
		tmp = lca(tochunkid(head(e)), tochunkid(tail(e)))
		if tolevel(tmp) == 1
			tmp = parent(tmp)
		end
		tmp == c.id
	end
	push!(c.added_edges, e)
	touch!(c)
end
@inline function delete_edge!(c::Chunk, e::AtomicEdge)
	@assert !(e in c.added_edges)
	@assert begin
		tmp = lca(tochunkid(head(e)), tochunkid(tail(e)))
		if tolevel(tmp) == 1
			tmp = parent(tmp)
		end
		tmp == c.id
	end
	push!(c.deleted_edges, e)
	touch!(c)
end

function Chunk(cgraph::ChunkedGraph, chunkid::ChunkID)
	m=MultiGraph{Label,AtomicEdge,CompositeEdgeSet}()
	return Chunk(cgraph, chunkid, Dict{Label, Vertex}(), m, Label(0))
end

@inline function tochunkid(c::Chunk)
	return c.id
end

@inline function tolevel(c::Chunk)
	return tolevel(c.id)
end

@inline function topos(c::Chunk)
	return topos(c.id)
end

@inline function isroot(c::Chunk)
	return isroot(c.id)
end

@inline function parent!(c::Chunk)
	return getchunk!(c.cgraph, parent(c.id))
end

function touch!(c::Void)
end
function touch!(c::Chunk)
	if c.clean
		c.clean = false
		c.modified = true
		touch!(c.parent)
	end
end

function uniquelabel!(c::Chunk)
	#TODO: Should erase clean flag?
	c.max_label += 1
	return tolabel(c.id, c.max_label)
end

#TODO: use a queue to maintain order of operations
function update!(c::Chunk)
	#TODO: tag updates with time since updates don't commute
	if c.clean
		return
	end

	# Update children first
	for child in c.children
		update!(child)
	end
	for child in c.children
		@assert child.clean
	end

	# Print debug messages
	# println("updating $(tolevel(c.id)), $(map(Int,topos(c.id))) V: +$(length(c.added_vertices))/-$(length(c.deleted_vertices)), E: +$(length(c.added_edges))/-$(length(c.deleted_edges))")
	global n_processed
	n_processed += 1
	# println("$n_processed/$(length(c.cgraph.chunks))")

	#vertices which need updates
	dirty_vertices = Set{Vertex}()
	redo_edge_sets = Set{CompositeEdgeSet}()

	# FIXME: We should upsize with the difference of added minus deleted
	upsize!(c.graph, length(c.added_vertices), length(c.added_edges))
	upsize!(c.vertices, length(c.added_vertices))

	# Insert added vertices
	# mark them as dirty_vertices as well
	for (l,v) in c.added_vertices
		@assert parent(tochunkid(v)) == c.id
		@assert v.parent == NULL_LABEL
		@assert !haskey(c.vertices, v.label)
		add_vertex!(c.graph, v.label)
		c.vertices[v.label] = v
		push!(dirty_vertices,v)
	end

	# Delete vertices and mark the vertices connected to 
	# the one we are deleting as dirty
	for (l,v) in c.deleted_vertices
		# for child in v.children
		# 	@assert get_vertex(c.cgraph, child).parent === NULL_LABEL
		# end

		for edge_set in incident_edges(c.graph, v.label)
			push!(redo_edge_sets, edge_set)
		end

		# this should delete all edges incident with v as well
		rem_vertex!(c.graph, v.label)
		delete!(c.vertices, v.label)
	end

	for e in redo_edge_sets
		for ee in revalidate!(c.cgraph, e)
			u = c.vertices[head(ee)]
			v = c.vertices[tail(ee)]
			add_edges!(c.graph,u.label,v.label,ee)
			push!(dirty_vertices, u, v)
		end
	end

	for child in c.children
		@assert child.clean
	end

	for ae in c.added_edges
		e=CompositeEdge(c.cgraph,ae)
		# TODO: Needs better handling
		u, v = c.vertices[head(e)], c.vertices[tail(e)]
		@assert tochunkid(u.label) != tochunkid(v.label) || (tolevel(u.label)==1 && tolevel(v.label)==1)
		@assert parent(tochunkid(u.label)) == parent(tochunkid(v.label)) == tochunkid(c)
		@assert haskey(c.vertices, u.label)
		@assert haskey(c.vertices, v.label)

		@assert haskey(c.graph.vertex_map, u.label)
		@assert haskey(c.graph.vertex_map, v.label)
		add_edges!(c.graph, u.label, v.label, e)
		push!(dirty_vertices, u, v)
	end

	for ae in c.deleted_edges
		e=CompositeEdge(c.cgraph,ae)
		u, v = c.vertices[head(e)], c.vertices[tail(e)]
		@assert isvalid(c.cgraph,e)
		@assert parent(tochunkid(head(e))) == parent(tochunkid(tail(e))) == tochunkid(c)
		@assert tochunkid(u.label) != tochunkid(v.label) || (tolevel(u.label)==1 && tolevel(v.label)==1)
		@assert parent(tochunkid(u.label)) == parent(tochunkid(v.label)) == tochunkid(c)
		rem_edges!(c.graph, u.label, v.label, e)
		push!(dirty_vertices, u, v)
	end

	if !isroot(c)
		#duplicated code because chain(dirty_vertices, c.deleted_vertices) is not type inferred
		#TODO: investigate or refactor
		for v in dirty_vertices
			if v.parent != NULL_LABEL
				@assert tochunkid(v.parent) == tochunkid(c)
				@assert parent(tochunkid(v.parent)) == tochunkid(c.parent)
				delete_vertex!(c.parent, v.parent)
				v.parent = NULL_LABEL
			end
		end
		for (l,v) in c.deleted_vertices
			if v.parent != NULL_LABEL
				@assert tochunkid(v.parent) == tochunkid(c)
				@assert parent(tochunkid(v.parent)) == tochunkid(c.parent)
				delete_vertex!(c.parent, v.parent)
				v.parent = NULL_LABEL
			end
		end
		c.parent.clean = false
	end

	cc = connected_components(c.graph, map(x->x.label, dirty_vertices))
	for component in cc
		if length(component) > 1
			l = uniquelabel!(c)
			new_vertex = Vertex(l, NULL_LABEL, component)
			for child_label in component
				c.vertices[child_label].parent = new_vertex.label
			end
			if !isroot(c)
				add_vertex!(c.parent, new_vertex)
			end
		else
			c.vertices[component[1]].parent = NULL_LABEL
		end
	end

	for v in values(c.vertices)
		@assert v.parent == NULL_LABEL || tochunkid(v.parent) == c.id
		@assert parent(tochunkid(v.label)) == c.id
	end



	empty!(c.added_edges)
	empty!(c.added_vertices)
	empty!(c.deleted_edges)
	empty!(c.deleted_vertices)
	c.clean = true
end
