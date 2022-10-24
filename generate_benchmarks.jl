using Benchmarks

# see doc for available options

function generate_benchmarks(args::Dict{String, Any})
    valid_std_systems = ["standard_bursting", "standard_lorenz", "bursting_limit_cycle", "lorenz_limit_cycle"]
    valid_ns_systems = ["ExplodingLorenz", "ShiftingLorenz", "ShrinkingLorenz", "PaperLorenzBigChange", "PaperLorenzSmallChange"]
    valid_trial_systems = ["trial_lorenz", "split_lorenz"]
    valid_regimes = ["bursting_neuron_regimes"]

    benchmark_system = args["name"]::String
    process_noise_level = args["process_noise_level"]::Float32
    num_T = args["num_T"]::Int
    ΔT = args["delta_T"]::Float32
    transient_T = args["transient_T"]::Int
    PLOT = args["PLOT"]::Bool
    plot_title = args["plot_title"]::String
    SAVE = args["SAVE"]::Bool
    save_dir = args["save_dir"]::String
    num_trials = args["num_trials"]::Int
    MARKERS = args["MARKERS"]::Bool
    p_change = args["p_change"]::Vector{String}
    lorenz_sys = args["lorenz_trial_sys"]::String

    if !any(x->benchmark_system in x, [valid_ns_systems,valid_trial_systems, valid_std_systems,valid_regimes])
        throw("$benchmark_system is not a valid system")
    end

    if benchmark_system in valid_std_systems
        std_3d_benchmark(benchmark_system, num_T=num_T, ΔT=ΔT, transient_T=transient_T, 
                        plot_title=plot_title, PLOT=PLOT, save_dir=save_dir, SAVE=SAVE,MARKER=MARKERS,
                        process_noise_level=process_noise_level)
    elseif benchmark_system in valid_regimes
        bursting_neuron_regimes(num_T=num_T, ΔT=ΔT, transient_T=transient_T,
                                PLOT=PLOT, save_dir=save_dir, SAVE=SAVE,
                                process_noise_level=process_noise_level)
    elseif benchmark_system in valid_ns_systems
        
        p_change_fun = Vector{Function}()
        for p_f in p_change
            p_sym = Symbol(p_f)
            p_fun = @eval $p_sym
            push!(p_change_fun, p_fun)
        end
        ns_3d_benchmark(benchmark_system, p_change=p_change_fun, num_T=num_T, ΔT=ΔT, transient_T=transient_T,
                        plot_title=plot_title,PLOT=PLOT, save_dir=save_dir,SAVE=SAVE, 
                        process_noise_level=process_noise_level)
    elseif benchmark_system in valid_trial_systems
        trial_benchmark("split_lorenz", num_trials, seq_length=num_T, ΔT=ΔT, transient_T=transient_T, 
                        plot_title=plot_title, PLOT=PLOT, save_dir=save_dir,SAVE=SAVE,
                        process_noise_level=process_noise_level ,lorenz_sys=lorenz_sys)
    else
        throw("something went wrong")
    end
end

generate_benchmarks(parse_commandline())