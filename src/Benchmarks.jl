module Benchmarks

include("ds_equations.jl")
export loop_burstn!, ns_bursting_neuron!, ns_lorenz!

include("ds_models.jl")
export lorenz, bursting_neuron

include("ns_models.jl")
export ns_benchmark_systems, ns_systems_bench

include("utils.jl")
export gen_path, save_series, generate_trajectories, sigmoid, linear, exponential

include("trial_ns.jl")
export trial_ns_lorenz,split_lorenz_trials

include("parsing.jl")
export parse_commandline, generate_benchmarks

include("data_gen.jl")
export std_3d_benchmark, bursting_neuron_regimes, ns_3d_benchmark, trial_benchmark


end
