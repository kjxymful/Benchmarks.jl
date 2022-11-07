using JSON
using ArgParse

function generate_benchmarks(args::Dict{String, Any})
    valid_std_systems = ["standard_bursting", "standard_lorenz", "bursting_limit_cycle", "lorenz_limit_cycle"]
    valid_ns_systems = ["RampUpBN","SuddenBurstBN","ExplodingLorenz", "ShiftingLorenz", "ShrinkingLorenz", "PaperLorenzBigChange", "PaperLorenzSmallChange"]
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
    p_change = args["p_change"]::String
    lorenz_sys = args["lorenz_trial_sys"]::String
    snapshots = args["snapshots"]

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
        

        p_sym = Symbol(p_change)
        p_change_fun = @eval $p_sym
        ns_3d_benchmark(benchmark_system, p_change=p_change_fun, num_T=num_T, ΔT=ΔT, transient_T=transient_T,
                        plot_title=plot_title,PLOT=PLOT, save_dir=save_dir,SAVE=SAVE, 
                        process_noise_level=process_noise_level, snapshots=snapshots)
    elseif benchmark_system in valid_trial_systems
        trial_benchmark("split_lorenz", num_trials, seq_length=num_T, ΔT=ΔT, transient_T=transient_T, 
                        plot_title=plot_title, PLOT=PLOT, save_dir=save_dir,SAVE=SAVE,
                        process_noise_level=process_noise_level ,lorenz_sys=lorenz_sys)
    else
        throw("something went wrong")
    end
end

"""
    parse_commandline()

Parses all commandline arguments for execution of `generate_benchmarks.jl`.
"""
function parse_commandline()
    settings = ArgParseSettings()
    defaults = load_defaults()

    @add_arg_table settings begin
        # meta
        "--name"
        help = "The benchmark name"
        arg_type = String
        default = defaults["name"] |> String

        "--process_noise_level"
        help = "The process noise level"
        arg_type = Float32
        default = defaults["process_noise_level"] |> Float32
        
        "--num_T"
        help = "number of timesteps"
        arg_type = Int
        default = defaults["num_T"] |> Int

        "--delta_T"
        help = "difference between timesteps"
        arg_type = Float32
        default = defaults["delta_T"] |>Float32
        
        "--transient_T"
        help = "timesteps discarded for startup"
        arg_type = Int
        default = defaults["transient_T"] |> Int

        "--PLOT"
        help = "whether a plot is done"
        arg_type = Bool
        default = defaults["PLOT"] |> Bool
        
        "--plot_title"
        help = "the title to give to the plot"
        arg_type = String
        default = defaults["plot_title"] |> String

        "--SAVE"
        help = "whether to save data"
        arg_type = Bool
        default = defaults["SAVE"] |> Bool

        "--save_dir"
        help = "where to save the data"
        arg_type = String
        default = defaults["save_dir"] |> String

        "--num_trials"
        help = "number of trials in case of trial benchmark"
        arg_type = Int
        default = defaults["num_trials"] |> Int

        "--MARKERS"
        help = "whether to show data points"
        arg_type = Bool
        default = defaults["MARKERS"] |> Bool

        "--p_change"
        help = "Type of change for ns parameters"
        arg_type = String
        default = defaults["p_change"] |>String

        "--lorenz_trial_sys"
        help = "the lorenz system used to split"
        arg_type = String
        default = defaults["lorenz_trial_sys"] |> String

        "--snapshots"
        help = "whether to plot snapshots"
        arg_type = Bool
        default = defaults["snapshots"] |> Bool
    end
    return parse_args(settings)
end

load_defaults() =
    convert_to_Float32(JSON.parsefile(joinpath(pwd(), "settings", "benchmark_defaults.json")))

function convert_to_Float32(dict::Dict)
    for (key, val) in dict
        dict[key] = val isa AbstractFloat ? Float32(val) : val
    end
    return dict
end

