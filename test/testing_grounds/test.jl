# include("../src/Benchmarks.jl")
# using .Benchmarks

using Benchmarks

function test_ns3d_stability(args::Dict{String,Any})
    num_test = 10

    benchmark_system = args["name"]::String
    process_noise_level = args["process_noise_level"]::Float32
    num_T = args["num_T"]::Int
    ΔT = args["delta_T"]::Float32
    transient_T = args["transient_T"]::Int
    PLOT = args["PLOT"]::Bool
    plot_title = args["plot_title"]::String
    SAVE = args["SAVE"]::Bool
    save_dir = args["save_dir"]::String
    p_change = args["p_change"]::Vector{String}

    SAVE = false

    p_change_fun = Vector{Function}()
    for p_f in p_change
        p_sym = Symbol(p_f)
        p_fun = @eval $p_sym
        push!(p_change_fun, p_fun)
    end

    u₀ = [0.5 for i in 1:3]
    for i in 1:num_test
        ns_3d_benchmark(benchmark_system, p_change=p_change_fun, num_T=num_T, ΔT=ΔT, transient_T=transient_T,
            plot_title=plot_title, PLOT=PLOT, save_dir=save_dir, SAVE=SAVE,
            process_noise_level=process_noise_level, u0=u₀,eval=true, eval_run=i)
        u₀.+= randn(3)
    end
end

test_ns3d_stability(parse_commandline())