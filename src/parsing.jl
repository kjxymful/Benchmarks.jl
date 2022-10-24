using JSON
using ArgParse

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
        arg_type = Vector{String}
        default = defaults["p_change"] |>Vector{String}

        "--lorenz_trial_sys"
        help = "the lorenz system used to split"
        arg_type = String
        default = defaults["lorenz_trial_sys"] |> String
    end
    return parse_args(settings)
end

load_defaults() =
    convert_to_Float32(JSON.parsefile(joinpath(pwd(), "settings", "defaults.json")))

function convert_to_Float32(dict::Dict)
    for (key, val) in dict
        dict[key] = val isa AbstractFloat ? Float32(val) : val
    end
    return dict
end