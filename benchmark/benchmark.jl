push!(LOAD_PATH, dirname(@__FILE__))
include("../src/chunkedgraphs/ChunkedGraphs.jl")
using ChunkedGraphs
using BenchmarkTools, Compat

settings = JSON.parsefile("benchmark.conf")

getvertexbench = @benchmarkable getvertex!(cgraph, 0x01040201000004ca) setup = begin
	cgraph = ChunkedGraph($(settings["graphpath"]))
end

getrootbench = @benchmarkable root!(cgraph, v) setup = begin
	cgraph = ChunkedGraph($(settings["graphpath"]))
	v = getvertex!(cgraph, 0x01040201000004ca)
end

getleavesbench = @benchmarkable leaves!(cgraph, p, 1) setup = begin
	cgraph = ChunkedGraph($(settings["graphpath"]))
	p = root!(cgraph, getvertex!(cgraph, 0x01040201000004ca))
end

run(getvertexbench)

run(getrootbench)

run(getleavesbench)
