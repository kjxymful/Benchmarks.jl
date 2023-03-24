module Benchmarks

using Statistics
using StatsBase
using DynamicalSystems

const valid_std_systems = ["standard_bursting", "standard_lorenz", "bursting_limit_cycle", "lorenz_limit_cycle"]
const valid_ns_systems = ["Roessler","EasyBN", "SimpleRampBN","RampUpBN", "StopBurstBN", "ExplodingLorenz", "ShiftingLorenz", "ShrinkingLorenz", "PaperLorenzBigChange", "PaperLorenzSmallChange"]
const valid_trial_systems = ["trial_lorenz", "split_lorenz"]
const valid_regimes = ["bursting_neuron_regimes"]

include("ds_equations.jl")

include("ds_models.jl")
export lorenz, bursting_neuron, ds_systems

include("ns_models.jl")
export ns_benchmark_systems, ns_systems_bench, ns_systems

include("utils.jl")
export generate_trajectories, standardize_to_nsseries,
        linear, sigmoid, exponential, complicated, complicated_lorenz,
        TP_loc, lorenz_at_t, BN_at_t, dynamical_noise_callback, rssl_at_t

include("trial_ns.jl")
export trial_ns_lorenz, split_lorenz_trials

include("data_gen.jl")

include("parsing.jl")
export parse_commandline, generate_benchmarks

end
