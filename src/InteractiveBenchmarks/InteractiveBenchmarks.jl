module InteractiveBenchmarks

using InteractiveDynamics, Makie, GLMakie, DynamicalSystems
using DifferentialEquations

export interactive_benchmarks
include("ds_pn_models.jl")
include("utillity.jl")
include("interactive.jl")

end