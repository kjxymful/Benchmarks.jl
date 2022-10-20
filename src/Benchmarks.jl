module Benchmarks

include("ds_models.jl")
export lorenz, bursting_neuron

include("ns_models.jl")
export ns_lorenz_systems

include("utils.jl")
export gen_path, save_series, gen_series, generate_ns_trajectories, sigmoid, linear, exponential

include("data_gen.jl")
export std_3d_benchmark, bursting_neuron_regimes, ns_3d_benchmark

end
