push!(LOAD_PATH, dirname(@__FILE__))
include("../../chunkedgraphs/ChunkedGraphs.jl")
using ChunkedGraphs
using CloudVolume

import PyCall
import JSON

settings = JSON.parsefile("graph.conf")

function buildgraphbase!(graph::ChunkedGraph, edgetasks::StorageWrapper, ranges::Tuple{StepRange,StepRange,StepRange})
	total = length(ranges[1]) * length(ranges[2]) * length(ranges[3])
	i = 0
	for z in ranges[3], y in ranges[2], x in ranges[1]
		i += 1
		chunkid = ChunkedGraphs.world_to_chunkid(graph, (x, y, z))
		prefix = "$(x)-$(x + step(ranges[1]))_$(y)-$(y + step(ranges[2]))_$(z)-$(z + step(ranges[3]))"
		atomicedges = edgetasks.val[:get_file]("$(prefix)_atomicedges.bin")
		atomicedges = Vector{AtomicEdge}(reinterpret(AtomicEdge, Vector{UInt8}(atomicedges)))

		rg2cg = edgetasks.val[:get_file]("$(prefix)_rg2cg.bin")
		rg2cg = Dict{UInt64, Label}(reinterpret(Pair{UInt64, Label}, Vector{UInt8}(rg2cg)))
		cg2rg = map(reverse, rg2cg)
		write("$(settings["graphpath"])/cg2rg_$(ChunkedGraphs.stringify(chunkid)).bin", collect(cg2rg))
		println("$i/$total | $prefix: Adding $(length(rg2cg)) vertices and $(length(atomicedges)) edges")

		add_atomic_vertices!(graph, collect(values(rg2cg)))
		add_atomic_edges!(graph, atomicedges)
	end

	update!(graph)
	save!(graph)
end

function buildgraphlayer!(graph::ChunkedGraph, edgetasks::StorageWrapper, ranges::Tuple{StepRange,StepRange,StepRange},
			interchunk_ranges::Tuple{StepRange,StepRange,StepRange})
	total = length(interchunk_ranges[1]) * length(ranges[2]) * length(ranges[3]) + 
			length(ranges[1]) * length(interchunk_ranges[2]) * length(ranges[3]) + 
			length(ranges[1]) * length(ranges[2]) * length(interchunk_ranges[3])

	@time begin
		i = 0
		for z in ranges[3], y in ranges[2], x in interchunk_ranges[1]
			i += 1
			prefix = "$(x)-$(x + 2)_$(y)-$(y + step(ranges[2]))_$(z)-$(z + step(ranges[3]))"
			atomicedges = edgetasks.val[:get_file]("$(prefix)_atomicedges.bin")
			atomicedges = Vector{AtomicEdge}(reinterpret(AtomicEdge, Vector{UInt8}(atomicedges)))
			println("$i/$total | $prefix: Adding $(length(atomicedges)) edges")
			add_atomic_edges!(graph, atomicedges)
		end
		for z in ranges[3], y in interchunk_ranges[2], x in ranges[1]
			i += 1
			prefix = "$(x)-$(x + step(ranges[1]))_$(y)-$(y + 2)_$(z)-$(z + step(ranges[3]))"
			atomicedges = edgetasks.val[:get_file]("$(prefix)_atomicedges.bin")
			atomicedges = Vector{AtomicEdge}(reinterpret(AtomicEdge, Vector{UInt8}(atomicedges)))
			println("$i/$total | $prefix: Adding $(length(atomicedges)) edges")
			add_atomic_edges!(graph, atomicedges)
		end
		for z in interchunk_ranges[3], y in ranges[2], x in ranges[1]
			i += 1
			prefix = "$(x)-$(x + step(ranges[1]))_$(y)-$(y + step(ranges[2]))_$(z)-$(z + 2)"
			atomicedges = edgetasks.val[:get_file]("$(prefix)_atomicedges.bin")
			atomicedges = Vector{AtomicEdge}(reinterpret(AtomicEdge, Vector{UInt8}(atomicedges)))
			println("$i/$total | $prefix: Adding $(length(atomicedges)) edges")
			add_atomic_edges!(graph, atomicedges)
		end

		update!(graph)
		save!(graph)
	end
end

function buildgraph(settings::Dict{String,Any})
	# Ensure graph directory exists
	if !isdir(settings["graphpath"])
		mkdir(settings["graphpath"])
	end

	# Init Cloud Storage and Chunked Graph
	graph = ChunkedGraph(settings["graphpath"], settings)
	edgetasks = StorageWrapper(settings["edgetasks"])

	ranges = ([settings["offset"][i] : settings["chunksize"][i] : settings["offset"][i] + settings["size"][i] - 1 for i in 1:3]...)

	# Build Base Layer
	buildgraphbase!(graph, edgetasks, ranges)

	newstep = settings["chunksize"]
	iterations = Int(ceil(log2(maximum(cld.(settings["offset"] .+ settings["size"], settings["chunksize"])))))
	for j = 1:iterations
		newstep = 2 * newstep
		newstart = @. newstep * cld(settings["offset"], newstep) + div(newstep, 2) - 1
		interchunk_ranges = ([newstart[i] : newstep[i] : settings["offset"][i] + settings["size"][i] for i in 1:3]...)
		buildgraphlayer!(graph, edgetasks, ranges, interchunk_ranges)
	end
end

function nextmult(n::Integer, m::Integer)
	return n >= 0 ? div((n + m - 1), m) * m : div(n, m) * m
end

function prevmult(n::Integer, m::Integer)
	return -nextmult(-n, m)
end

roisize = prevpow2.(last.(ranges))

minpoint = (settings["offset"]...)
midpoint = prevpow2.(last.(ranges))
maxpoint = ((settings["offset"] .+ settings["size"])...)


function buildgraph_recursive(settings::Dict{String, Any}, bbox::Tuple{StepRange{Int64},StepRange{Int64},StepRange{Int64}})
	if any(last.(bbox) .< (settings["offset"]...)) ||
			any(start.(bbox)) .> ((settings["offset"] .+ settings["size"])...)
		return
	end



end

# Dataset origin and World origin are usually not the same. We keep it that way
# as it makes converting chunks to coordinates a bit easier...
function buildgraph(settings::Dict{String, Any})
	chunksize = (settings["step"]...)
	datasetsize = (settings["size"]...)
	datasetbbox = range.((settings["offset"]...), datasetsize)
	worldbbox = range.((0, 0, 0), nextpow2.(last.(datasetbbox)))
	depth = Int(log2(maximum(div.(last.(worldbbox), chunksize))))

	halfsize = 2^(depth-1) .* chunksize
	offsets = first.(worldbbox)

	buildgraph_recursive(range.(offsets, halfsize, 2))



	midpoint = div.(last.(ranges) .- first.(ranges), 2)
	
	right = nextpow2.(last.(ranges))
	mid = prevpow2.(r .- 1)
	left = r .- 2 .* m
end