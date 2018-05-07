import Base: show

show(io::IO, cgraph::ChunkedGraph) = print("ChunkedGraph with $(length(cgraph.chunks)) chunks in memory")

function loadchunk(cgraph::ChunkedGraph, chunkid::ChunkID)
	vertex_map = Dict{Label, Vertex}()
	mgraph = MultiGraph()
	max_label = Label(0)

	prefix = stringify(chunkid)
	path = expanduser(joinpath(cgraph.path, String("$(prefix).chunk")))
	if !isfile(path)
		return Chunk(cgraph, chunkid, vertex_map, mgraph, max_label)
	end

	print("Loading from $path...")
	f = open(path, "r")

	# Check File Version
	version = read(f, UInt64)
	@assert version === UInt64(1)

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
		(u, v, atomic_edge_cnt) = read(f, UInt64, 3)
		atomic_edges = read(f, AtomicEdge, atomic_edge_cnt)
		add_edges!(mgraph, u, v, atomic_edges)
	end

	close(f)
	println("done.")
	return Chunk(cgraph, chunkid, vertex_map, mgraph, max_label)
end

function getchunk!(cgraph::ChunkedGraph, chunkid::ChunkID)
	if !haskey(cgraph.chunks, chunkid)
		cgraph.chunks[chunkid] = loadchunk(cgraph, chunkid)
		cgraph.lastused[chunkid] = 0
	end

	c = cgraph.chunks[chunkid]::Chunk
	if cgraph.evictionmode === false
		gentletouch!(c)
	end
	return c
end

function getvertex!(cgraph::ChunkedGraph, l::Label)
	return getchunk!(cgraph, parent(cgraph, tochunkid(l))).vertices[l]
end

function hasvertex!(cgraph::ChunkedGraph, l::Label)
	return haskey(getchunk!(cgraph, parent(cgraph, tochunkid(l))).vertices, l)
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
	print("Saving to $(path)...")

	buf = IOBuffer()
	write(buf, UInt64(1)) # Version
	write(buf, UInt64(c.max_label)) # Max SegID
	write(buf, UInt64(length(c.vertices))) # Vertex Count
	write(buf, UInt64(length(LightGraphs.edges(c.graph.graph)))) # Edge Count

	for vertex in values(c.vertices)
		write(buf, UInt64(vertex.label)) # Vertex Label
		write(buf, UInt64(vertex.parent)) # Vertex Parent
		write(buf, UInt64(length(vertex.children))) # Vertex Children Count
		write(buf, convert(Vector{UInt64}, vertex.children))
	end

	for (edge, atomic_edges) in c.graph.edge_map
		write(buf, UInt64(c.graph.inverse_vertex_map[edge[1]]))
		write(buf, UInt64(c.graph.inverse_vertex_map[edge[2]]))
		write(buf, UInt64(length(atomic_edges)))
		write(buf, convert(Vector{AtomicEdge}, collect(atomic_edges)))
	end

	f = open(path, "w")
	write(f, buf.data)
	close(f)
	close(buf)
	println("done.")
	
	c.modified = false
end

function run_eviction!(cgraph::ChunkedGraph)
	cgraph.evictionmode = true
	while length(cgraph.chunks) > div(cgraph.CACHESIZE, 2)
		convict, priority = DataStructures.peek(cgraph.lastused)
		children = cgraph.chunks[convict].children
		if length(children) > 0
			# Spare him, kill the children first
			# He'll be put back on death row when his children die
			dequeue!(cgraph.lastused, convict)
			continue
		else
			c = cgraph.chunks[convict]
			@assert !isroot(c)
			@assert c.id != c.cgraph.SECOND_ID

			# Make sure all queued vertices of parent chunk have been generated.
			update!(c.parent)
			evict!(cgraph.chunks[convict])
		end
	end
	cgraph.evictionmode = false
end

function evict!(c::Chunk)
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

	if length(c.parent.children) === 0 && !haskey(c.cgraph.lastused, c.parent.id)
		c.cgraph.lastused[c.parent.id] = priority + 1
	end
end
