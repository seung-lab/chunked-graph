import Base: show, write, read

# using Logging
# const iologger = Logger("iologger", level = OFF)

show(io::IO, cgraph::ChunkedGraph) = print("ChunkedGraph with $(length(cgraph.chunks)) chunks in memory")

function write(io::IO, e::CompositeEdgeSet)
	write(io, e.u)
	write(io, e.v)
	write(io, e.max_affinity)
	write(io, e.nonempty)
	write(io, UInt64(length(e.children)))
	for c in e.children
		write(io, c)
	end
end

function read(io::IO, ::Type{CompositeEdgeSet})
	u = read(io, Label)
	v = read(io, Label)
	aff = read(io, Affinity)
	nonempty = read(io, Bool)
	n = read(io, UInt64)
	children = CompositeEdgeSet[read(io, CompositeEdgeSet) for i in 1:n]
	return CompositeEdgeSet(u, v, aff, children, nonempty)
end

function write(io::IO, e::AtomicEdge)
	write(io, e.u)
	write(io, e.v)
	write(io, e.affinity)
end

function read(io::IO, ::Type{AtomicEdge})
	return AtomicEdge(read(io, Label), read(io, Label), read(io, Affinity))
end

function loadchunk(cgraph::ChunkedGraph, chunkid::ChunkID)
	vertex_map = Dict{Label, Vertex}()
	mgraph = MultiGraph{Label, AtomicEdge, CompositeEdgeSet}()
	max_label = Label(0)

	prefix = stringify(chunkid)
	path = expanduser(joinpath(cgraph.path, String("$(prefix).chunk")))
	if !isfile(path)
		return Chunk(cgraph, chunkid, vertex_map, mgraph, max_label)
	end
	
	f = open(path, "r")

	# Check File Version
	version = read(f, UInt64)
	@assert version === UInt64(2)
	
	# Read Chunk Info
	(max_label, v_cnt, e_cnt) = read(f, UInt64, 3)

	# Allocate Graph
	sizehint!(mgraph, floor(UInt32, 1.5 * v_cnt), floor(UInt32, 1.5 * e_cnt))
	sizehint!(vertex_map, floor(UInt32, 1.5 * v_cnt))
	
	# Read Vertices
	for i in range(1, v_cnt)
		(label, parent, child_cnt) = read(f, UInt64, 3)
		children = read(f, UInt64, child_cnt)
		add_vertex!(mgraph, label)
		vertex_map[label] = Vertex(label, parent, children)
	end
	# Read EdgeMap
	for i in range(1, e_cnt)
		edge_set = read(f, CompositeEdgeSet)
		add_edges!(mgraph, head(edge_set), tail(edge_set), edge_set)
	end

	close(f)
	return Chunk(cgraph, chunkid, vertex_map, mgraph, max_label)
end

function getindex!(cgraph::ChunkedGraph, chunkid::ChunkID)
	return getchunk!(cgraph, chunkid)
end

@inline function gentle_touch!(c::Chunk)
	c.cgraph.lastused[c.id] = time_ns()
end

function getchunk!(cgraph::ChunkedGraph, chunkid::ChunkID)
	if !haskey(cgraph.chunks, chunkid)
		tmp = loadchunk(cgraph, chunkid)
		if tmp == nothing
			tmp = Chunk(cgraph, chunkid)
		end
		cgraph.chunks[chunkid] = tmp::Chunk
		cgraph.lastused[chunkid] = 0
	end

	c = cgraph.chunks[chunkid]::Chunk
	if !cgraph.eviction_mode
		gentle_touch!(c::Chunk)
	end
	return c
end

function run_eviction!(cgraph::ChunkedGraph)
	cgraph.eviction_mode = true
	while length(cgraph.chunks) > CACHESIZE
		convict, priority = DataStructures.peek(cgraph.lastused)
		children = cgraph.chunks[convict].children
		if length(children) > 0
			#spare him, kill the children first
			#he'll be put back on death row when his children die
			dequeue!(cgraph.lastused, convict)
			continue
		else
			evict!(cgraph.chunks[convict])
		end
	end
	cgraph.eviction_mode = false
end

function getvertex!(cgraph::ChunkedGraph, l::Label)
	return getchunk!(cgraph, parent(tochunkid(l))).vertices[l]
end

function hasvertex!(cgraph::ChunkedGraph, l::Label)
	return haskey(getchunk!(cgraph, parent(tochunkid(l))).vertices,l)
end

function save!(cgraph::ChunkedGraph, force::Bool = false)
	for c in collect(values(cgraph.chunks))
		if c.modified || force
			save!(c)
		end
	end
end

function save!(c::Chunk)
	@assert c.clean

	prefix = stringify(c.id)
	path = expanduser(joinpath(c.cgraph.path, "$(prefix).chunk"))
	# info(iologger, "Saving to $(path)...")
	
	buf = IOBuffer()
	write(buf, UInt64(2)) # Version
	write(buf, UInt64(c.max_label)) # Max SegID
	write(buf, UInt64(length(c.vertices))) # Vertex Count
	write(buf, UInt64(length(LightGraphs.edges(c.graph.graph)))) # Edge Count
	for vertex in values(c.vertices)
		write(buf, UInt64(vertex.label)) # Vertex Label
		write(buf, UInt64(vertex.parent)) # Vertex Parent
		write(buf, UInt64(length(vertex.children))) # Vertex Children Count
		write(buf, convert(Vector{UInt64}, vertex.children))
	end
		
	for edgeset in values(c.graph.edge_map)
		write(buf, edgeset)
	end

	f = open(path, "w")
	write(f, buf.data)
	close(f)
	close(buf)
	# info(iologger, "Done saving to $(path)")
	
	c.modified = false
end

const eviction_count = DataStructures.DefaultDict{ChunkID,Int}(0)
function evict!(c::Chunk)
	global eviction_count
	t = (eviction_count[c.id] += 1)
	# if t > 1
	# 	warn(iologger, "evicted $(stringify(c.id)) $(t) times")
	# end

	@assert !isroot(c)
	@assert c.id != SECOND_ID
	update!(c)
	while !isempty(c.children)
		evict!(c.children[end])
	end

	if c.modified
		save!(c)
	end

	filter!(x->x != c, c.parent.children)
	delete!(c.cgraph.chunks, c.id)
	priority = c.cgraph.lastused[c.id]
	dequeue!(c.cgraph.lastused, c.id)
	if length(c.parent.children) == 0 && !haskey(c.cgraph.lastused, c.parent.id)
		c.cgraph.lastused[c.parent.id] = priority + 1
	end
	unlock(c.flock)
	close(c.flock)
end
