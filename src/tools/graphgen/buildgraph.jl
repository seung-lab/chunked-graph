addprocs(32)
@everywhere include("../../chunkedgraphs/ChunkedGraphs.jl")
@everywhere using ChunkedGraphs
@everywhere using CloudVolume
@everywhere using Retry

import JSON

settings = JSON.parsefile("pinky40.conf")

@everywhere function buildgraphbase!(settings::Dict{String,Any}, graph::ChunkedGraph, edgetasks::StorageWrapper, ranges::Tuple{StepRange,StepRange,StepRange})
	total = length(ranges[1]) * length(ranges[2]) * length(ranges[3])
	i = 0
	for z in ranges[3], y in ranges[2], x in ranges[1]
		i += 1
		chunkid = ChunkedGraphs.world_to_chunkid(graph, (x, y, z))
		prefix = "$(x)-$(x + step(ranges[1]))_$(y)-$(y + step(ranges[2]))_$(z)-$(z + step(ranges[3]))"
		atomicedges = nothing
		@repeat 7 try
			atomicedges = edgetasks.val[:get_file]("$(prefix)_atomicedges.bin")
			atomicedges = Vector{AtomicEdge}(reinterpret(AtomicEdge, Vector{UInt8}(atomicedges)))
		catch e
			println(e)
			println("Error while processing $(prefix)_atomicedges.bin")
			@delay_retry if true === true end
		end

		rg2cg = nothing
		@repeat 7 try
			rg2cg = edgetasks.val[:get_file]("$(prefix)_rg2cg.bin")
			rg2cg = Dict{UInt64, Label}(reinterpret(Pair{UInt64, Label}, Vector{UInt8}(rg2cg)))
		catch
			println(e)
			println("Error while processing $(prefix)_rg2cg.bin")
			@delay_retry if true === true end
		end
		cg2rg = map(reverse, rg2cg)
		write("$(settings["graphpath"])/cg2rg_$(ChunkedGraphs.stringify(chunkid)).bin", collect(cg2rg))
		println("$i/$total | $prefix: Adding $(length(rg2cg)) vertices and $(length(atomicedges)) edges")

		add_atomic_vertices!(graph, collect(values(rg2cg)))
		add_atomic_edges!(graph, atomicedges)
	end

	update!(graph)
	save!(graph)
end

@everywhere function buildgraphlayer!(graph::ChunkedGraph, edgetasks::StorageWrapper, ranges::Tuple{StepRange,StepRange,StepRange},
			interchunk_ranges::Tuple{StepRange,StepRange,StepRange})
	total = length(interchunk_ranges[1]) * length(ranges[2]) * length(ranges[3]) + 
			length(ranges[1]) * length(interchunk_ranges[2]) * length(ranges[3]) + 
			length(ranges[1]) * length(ranges[2]) * length(interchunk_ranges[3])

	println("Building Subgraph using ranges: $ranges and $interchunk_ranges")
	@time begin
		i = 0
		for z in ranges[3], y in ranges[2], x in interchunk_ranges[1]
			i += 1
			prefix = "$(x)-$(x + 2)_$(y)-$(y + step(ranges[2]))_$(z)-$(z + step(ranges[3]))"
			atomicedges = nothing
			@repeat 7 try
				atomicedges = edgetasks.val[:get_file]("$(prefix)_atomicedges.bin")
				atomicedges = Vector{AtomicEdge}(reinterpret(AtomicEdge, Vector{UInt8}(atomicedges)))
			catch e
				println(e)
				println("Error while processing $(prefix)_atomicedges.bin")
				@delay_retry if true === true end
			end
			println("$i/$total | $prefix: Adding $(length(atomicedges)) edges")
			add_atomic_edges!(graph, atomicedges)
		end
		for z in ranges[3], y in interchunk_ranges[2], x in ranges[1]
			i += 1
			prefix = "$(x)-$(x + step(ranges[1]))_$(y)-$(y + 2)_$(z)-$(z + step(ranges[3]))"
			atomicedges = nothing
			@repeat 7 try
				atomicedges = edgetasks.val[:get_file]("$(prefix)_atomicedges.bin")
				atomicedges = Vector{AtomicEdge}(reinterpret(AtomicEdge, Vector{UInt8}(atomicedges)))
			catch e
				println(e)
				println("Error while processing $(prefix)_atomicedges.bin")
				@delay_retry if true === true end
			end
			println("$i/$total | $prefix: Adding $(length(atomicedges)) edges")
			add_atomic_edges!(graph, atomicedges)
		end
		for z in interchunk_ranges[3], y in ranges[2], x in ranges[1]
			i += 1
			prefix = "$(x)-$(x + step(ranges[1]))_$(y)-$(y + step(ranges[2]))_$(z)-$(z + 2)"
			atomicedges = nothing
			@repeat 7 try
				atomicedges = edgetasks.val[:get_file]("$(prefix)_atomicedges.bin")
				atomicedges = Vector{AtomicEdge}(reinterpret(AtomicEdge, Vector{UInt8}(atomicedges)))
			catch e
				println(e)
				println("Error while processing $(prefix)_atomicedges.bin")
				@delay_retry if true === true end
			end
			println("$i/$total | $prefix: Adding $(length(atomicedges)) edges")
			add_atomic_edges!(graph, atomicedges)
		end

		update!(graph)
		save!(graph)
	end
end

@everywhere function buildgraph(settings::Dict{String,Any})
	# Ensure graph directory exists
	if !isdir(settings["graphpath"])
		mkdir(settings["graphpath"])
	end

	# Init Cloud Storage and Chunked Graph
	graph = ChunkedGraph(settings["graphpath"], settings)
	edgetasks = StorageWrapper(settings["edgetasks"])

	ranges = ([settings["offset"][i] : settings["chunksize"][i] : settings["offset"][i] + settings["size"][i] - 1 for i in 1:3]...)

	# Build Base Layer
	buildgraphbase!(settings, graph, edgetasks, ranges)

	newstep = settings["chunksize"]
	iterations = 3 #Int(ceil(log2(maximum(cld.(settings["size"], settings["chunksize"]))))) - 1
	for j = 1:iterations
		newstart = @. newstep + 2 * newstep * div((settings["offset"] + newstep), 2 * newstep) - 1
		interchunk_ranges = ([newstart[i] : 2 * newstep[i] : settings["offset"][i] + settings["size"][i] - 2 for i in 1:3]...)
		newstep = 2 * newstep
		buildgraphlayer!(graph, edgetasks, ranges, interchunk_ranges)
	end
end

@everywhere function nextmult(n::Integer, m::Integer)
	return n >= 0 ? div((n + m - 1), m) * m : div(n, m) * m
end

@everywhere function prevmult(n::Integer, m::Integer)
	return -nextmult(-n, m)
end

function buildgraph(settings::Dict{String,Any}, subdivisions::Integer=2)
	ranges = ([settings["offset"][i] : settings["chunksize"][i] : settings["offset"][i] + settings["size"][i] - 1 for i in 1:3]...)

	global_size = (nextpow2.(settings["size"])...)
	global_max = ((settings["offset"] .+ settings["size"])...)
	global_max_pow2 = nextpow2.(global_max)
	global_min = global_max_pow2 .- global_size
	block_count = 2^subdivisions
	block_size = maximum(div.(global_size, settings["chunksize"])) * div.(settings["chunksize"], block_count)

	new_ranges = collect.(range.(global_min, block_size, block_count))
	new_settings = Vector{Dict{String,Any}}()

	for x in new_ranges[1]
		for y in new_ranges[2]
			for z in new_ranges[3]
				s = deepcopy(settings)
				s["offset"] = max.(settings["offset"], [x, y, z])
				if any(s["offset"] .>= global_max)
					continue
				end
				block_idx = div.(s["offset"], block_size)
				s["size"] = min.(nextmult.(s["offset"] .+ 1, block_size) .- s["offset"], global_max .- s["offset"]) .- 1 #Don't connect last edge!
				s["graphpath"] = "$(settings["graphpath"])/$(subdivisions)_$(block_idx[1])_$(block_idx[2])_$(block_idx[3])"
				s["maxdepth"] = settings["maxdepth"] #- subdivisions + 1
				push!(new_settings, s)
			end
		end
	end

	pmap(buildgraph, new_settings)

	# TODO: Move to common folder

	# TODO: Modify graph.conf
	
	# TODO: Run remaining edges

end
#buildgraph(settings, 3)
