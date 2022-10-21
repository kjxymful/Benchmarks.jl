using Benchmarks

# see doc for available options

# # create standard benchmarks
# std_3d_benchmark("standard_lorenz")

# create data for different bursting neuron regimes
bursting_neuron_regimes(process_noise_level=0.6)

# # create non stationary benchmark systems
# # so far only lorenz systems were implemented
# ns_3d_benchmark("PaperLorenzSmallChange", process_noise_level=0.01)
