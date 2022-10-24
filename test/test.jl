include("../src/Benchmarks.jl")
using .Benchmarks

s = :linear
st = ["linear","linear"]
sym = Symbol.(st)

a = []
for l in st
    l_sym = Symbol(l)
    l_f = @eval $l_sym
    @show l_f
    push!(a,l_f)
end
@show a