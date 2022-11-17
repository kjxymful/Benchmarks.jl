module Benchmarks

using DynamicalSystems
using Statistics
using StatsBase

using Reexport

const valid_std_systems = ["standard_bursting", "standard_lorenz", "bursting_limit_cycle", "lorenz_limit_cycle"]
const valid_ns_systems = ["RampUpBN", "StopBurstBN", "ExplodingLorenz", "ShiftingLorenz", "ShrinkingLorenz", "PaperLorenzBigChange", "PaperLorenzSmallChange"]
const valid_trial_systems = ["trial_lorenz", "split_lorenz"]
const valid_regimes = ["bursting_neuron_regimes"]

include("ds_equations.jl")

include("ds_models.jl")
export lorenz, bursting_neuron, ds_systems

include("ns_models.jl")
export ns_benchmark_systems, ns_systems_bench, ns_systems

include("utils.jl")
export generate_trajectories

include("trial_ns.jl")
export trial_ns_lorenz, split_lorenz_trials

include("data_gen.jl")

include("parsing.jl")
export parse_commandline, generate_benchmarks

include("InteractiveBenchmarks/InteractiveBenchmarks.jl")
@reexport using .InteractiveBenchmarks

end
