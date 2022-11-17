using Benchmarks
using Test

valid_std_systems = ["standard_bursting", "standard_lorenz", "bursting_limit_cycle", "lorenz_limit_cycle"]
valid_ns_systems = ["RampUpBN", "StopBurstBN", "ExplodingLorenz", "ShiftingLorenz", "ShrinkingLorenz", "PaperLorenzBigChange", "PaperLorenzSmallChange"]
valid_trial_systems = ["trial_lorenz", "split_lorenz"]
valid_regimes = ["bursting_neuron_regimes"]

function wait_for_key(; prompt="press any key", io=stdin)
    setraw!(raw) = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid}, Int32), io.handle, raw)
    print(io, prompt)
    setraw!(true)
    read(io, 1)
    setraw!(false)
    nothing
end

@testset "Base Benchmarks" begin
    args = parse_commandline()

    all_sys = valid_std_systems
    append!(all_sys, valid_ns_systems)
    append!(all_sys, valid_regimes)

    T = args["T"]
    transient_T = args["transient_T"]
    for sys in all_sys
        args["name"] = sys
        @show sys
        if occursin("BN", sys) || occursin("bursting", sys)
            args["T"] = T * 10
            args["transient_T"] = transient_T * 10
        else
            args["T"] = T
            args["transient_T"] = transient_T
        end
        generate_benchmarks(args)
    end
    wait_for_key()
    rm("Figures/")
    rm("data/")
end

@testset "Dynamical Noise" begin
    args = parse_commandline()

    args["process_noise_level"] = 0.01
    all_sys = valid_std_systems
    append!(all_sys, valid_ns_systems)
    append!(all_sys, valid_regimes)

    T = args["T"]
    transient_T = args["transient_T"]
    for sys in all_sys
        args["name"] = sys
        @show sys
        if occursin("BN", sys) || occursin("bursting", sys)
            args["T"] = T * 10
            args["transient_T"] = transient_T * 10
        else
            args["T"] = T
            args["transient_T"] = transient_T
        end
        generate_benchmarks(args)
    end
    wait_for_key()
    rm("Figures/")
    rm("data/")
end
@testset "Trials" begin
    args = parse_commandline()

    all_sys = valid_std_systems
    append!(all_sys, valid_ns_systems)
    append!(all_sys, valid_regimes)

    T = args["T"]
    transient_T = args["transient_T"]
    for sys in valid_trial_systems
        args["name"] = sys
        generate_benchmarks(args)
    end
    rm("Figures/")
    rm("data/")
end
