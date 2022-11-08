using Benchmarks
using Test

function wait_for_key(; prompt="press any key", io=stdin)
    setraw!(raw) = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid}, Int32), io.handle, raw)
    print(io, prompt)
    setraw!(true)
    read(io, 1)
    setraw!(false)
    nothing
end

@testset "Benchmarks.jl" begin
    args = parse_commandline()
    valid_std_systems = ["standard_bursting", "standard_lorenz", "bursting_limit_cycle", "lorenz_limit_cycle"]
    valid_ns_systems = ["RampUpBN", "SuddenBurstBN", "ExplodingLorenz", "ShiftingLorenz", "ShrinkingLorenz", "PaperLorenzBigChange", "PaperLorenzSmallChange"]
    valid_trial_systems = ["trial_lorenz", "split_lorenz"]
    valid_regimes = ["bursting_neuron_regimes"]

    all_sys = valid_std_systems
    append!(all_sys, valid_ns_systems)
    append!(all_sys, valid_regimes)

    num_T = args["num_T"]
    transient_T = args["transient_T"]
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
    rm("data/", recursive=true)
    rm("Figures/", recursive=true)

    @testset "Trials" begin
        for sys in valid_trial_systems
            args["name"] = sys
            generate_benchmarks(args)
            rm("data/", recursive=true)
            rm("Figures/", recursive=true)
        end
    end
end
