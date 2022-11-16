using JSON
using ArgParse


function generate_benchmarks(args::Dict{String,Any})

    benchmark_system = args["name"]::String
    process_noise_level = args["process_noise_level"]::AbstractFloat
    num_T = args["num_T"]::Int
    Δt = args["delta_T"]::AbstractFloat
    transient_T = args["transient_T"]::Int
    PLOT = args["PLOT"]::Bool
    plot_title = args["plot_title"]::String
    SAVE = args["SAVE"]::Bool
    save_dir = args["save_dir"]::String
    num_trials = args["num_trials"]::Int
    MARKER = args["MARKERS"]::Bool
    p_change = args["p_change"]::String
    lorenz_sys = args["lorenz_trial_sys"]::String
    snapshots = args["snapshots"]
    plot_params = args["plot_params"]
    plot_name = args["exp_name"]

    check_validity(benchmark_system)
    if benchmark_system in valid_std_systems
        std_3d_benchmark(benchmark_system; num_T, Δt,transient_T,plot_title,PLOT, save_dir, SAVE, MARKER,
            process_noise_level, plot_name)
    elseif benchmark_system in valid_regimes
        bursting_neuron_regimes(;num_T, Δt,transient_T,PLOT, save_dir,SAVE,
            process_noise_level, plot_name)
    elseif benchmark_system in valid_ns_systems
        p_sym = Symbol(p_change)
        p_change = @eval $p_sym
        ns_3d_benchmark(benchmark_system; p_change, num_T, Δt,transient_T,
            plot_title, PLOT, save_dir, SAVE,
            process_noise_level, snapshots,
            plot_params, plot_name)
    elseif benchmark_system in valid_trial_systems
        trial_benchmark("split_lorenz", num_trials; seq_length=num_T, Δt, transient_T,
            plot_title, PLOT, save_dir, SAVE,
            process_noise_level, lorenz_sys)
    else
        throw("something went wrong")
    end
end



"""
    parse_commandline()

Parses all commandline arguments for execution of `generate_benchmarks.jl`.
"""
function parse_commandline(;path="")
    settings = ArgParseSettings()
    defaults = load_defaults(path)

    @add_arg_table! settings begin
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

        "--plot_params"
        help = "whether to plot parameter change"
        arg_type = Bool
        default = defaults["plot_params"] |> Bool

        "--exp_name"
        help = "name as which to save the plot"
        arg_type = String
        default = defaults["exp_name"] |> String
        
        "--interactive"
        help="starts the interactive exploration framework"
        arg_type = Bool
        default = defaults["interactive"] |> Bool

        "--interactive_noise"
        help="adds process noise as slider"
        arg_type = Float32
        default = defaults["interactive_noise"] |> Float32
    end
    return parse_args(settings)
end

load_defaults(path::String) =
    convert_to_Float32(JSON.parsefile(joinpath(pwd(), path*"settings", "benchmark_defaults.json")))

function convert_to_Float32(dict::Dict)
    for (key, val) in dict
        dict[key] = val isa AbstractFloat ? Float32(val) : val
    end
    return dict
end

