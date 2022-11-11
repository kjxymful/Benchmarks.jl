using Benchmarks
using DynamicalSystems

sys, params = ns_benchmark_systems("ShrinkingLorenz", linear, 150.0f0, zeros(3), transient_T=50)
ts = trajectory(sys, 150, Ttr=50)
p = plot3d(ts[:, 1], ts[:, 2], ts[:, 3], border=:none, axis=nothing, legend=nothing, lc=(:viridis), linez=collect(0:0.01:150))
savefig("examples/head.png")


args = parse_commandline()
valid_std_systems = ["standard_bursting", "standard_lorenz", "bursting_limit_cycle", "lorenz_limit_cycle"]
valid_ns_systems = ["RampUpBN", "StopBurstBN", "ExplodingLorenz", "ShiftingLorenz", "ShrinkingLorenz", "PaperLorenzBigChange", "PaperLorenzSmallChange"]
valid_trial_systems = ["trial_lorenz", "split_lorenz"]
valid_regimes = ["bursting_neuron_regimes"]

all_sys = valid_std_systems
append!(all_sys, valid_ns_systems)
append!(all_sys, valid_regimes)

num_T = args["num_T"]
transient_T = args["transient_T"]
args["SAVE"] = false
for sys in all_sys
    args["name"] = sys
    if occursin("BN", sys) || occursin("bursting", sys)
        args["num_T"] = num_T * 10
        args["transient_T"] = transient_T * 10
    else
        args["num_T"] = num_T
        args["transient_T"] = transient_T
    end
    generate_benchmarks(args)
end


for sys in valid_trial_systems
    args["name"] = sys
    generate_benchmarks(args)
end