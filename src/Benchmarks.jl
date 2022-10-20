module Benchmarks

include("ds_models.jl")
export lorenz, bursting_neuron

include("utils.jl")
export gen_path, save_series, gen_series

include("data_gen.jl")
export std_3d_benchmark, bursting_neuron_regimes

end
