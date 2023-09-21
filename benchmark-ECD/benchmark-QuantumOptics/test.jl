

include("benchmarkutils.jl")
testdic = Dict("N"=>2, "t"=>3)

@__DIR__
@__FILE__
benchmarkutils.save("blub", testdic)
const rootpath = abspath(".")