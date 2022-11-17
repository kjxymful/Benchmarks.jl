module InteractiveBenchmarks

using InteractiveDynamics, GLMakie, DynamicalSystems
using OrdinaryDiffEq

export interactive_benchmarks
include("ds_pn_models.jl")
include("utillity.jl")
include("interactive.jl")

end