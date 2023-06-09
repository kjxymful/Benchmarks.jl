using Plots
using NPZ: npzwrite

"""
    std_3d_benchmark(System; T=150, Δt=0.01, transient_T=20, process_noise_level=0.0, PLOT=true, plot_title="", plot_name="", MARKER=false, save_dir="", SAVE=true)

create a time series of the specified Benchmark System

implemented Benchmark Systems:
- standard_bursting: Bursting Neuron as used in all Durstewitz Lab papers
- standard_lorenz: Lorenz, with ρ = 28
- lorenz_limit_cycle: A limit cycle of the Lorenz System with ρ=24
- bursting_limit_cycle: The Bursting neuron limit cycle at gₙₘ₀ₐ=10

Parameter
---------
System : One of the valid Benchmark Systems

Kwargs
------
num_T : number of timesteps, default 15000 (Int)

ΔT : time between steps, default 0.01 (Float)

transient_T : timesteps used as start up to avoid transients, default 2000 (Int)

plot_title : Title for the plot, default: same as System (String)

PLOT : Specifies whether a plot is done, default true (BOOL)

save_dir : directory to save the data in, default data/benchmarks (String)

SAVE : Specifies whether the timeseries is saved or not (BOOL)

MARKER: Show data points in plot (BOOL)

process_noise_level : The ratio of process noise, default 0 (Float)

plot_name : the name of the file as which plot and data are saved

Returns
-------
time series : Time series of the specified system (AbstractMatrix)

"""
function std_3d_benchmark(System::String; T=150, Δt=0.01, transient_T=20, process_noise_level=0.0, PLOT=true, plot_title="", plot_name="", MARKER=false, save_dir="", SAVE=true)
    @assert System in valid_std_systems "$System is not a valid System: $valid_std_systems"
    exp_name = isempty(plot_name) ? String(System) : plot_name

    if System == "standard_bursting"
        ds::ContinuousDynamicalSystem = bursting_neuron()
        μ = 10.2
    elseif System == "standard_lorenz"
        μ = 28.0
        ds = lorenz()
    elseif System == "bursting_limit_cycle"
        ds = bursting_neuron(u0=[-60.0, 0.0386, 0.0231], gₙₘ₀ₐ=10.0)
        μ = 10.1
    elseif System == "lorenz_limit_cycle"
        μ = 23.0
        ds = lorenz(ρ=23.0)
    end
    tseries = generate_trajectories(ds, T, transient_T;Δt, process_noise_level, PLOT=false, model_name=System)

    time = 0:Δt:T

    if PLOT
        marker = MARKER ? :d : false
        p = plot3d(tseries[:, 1], tseries[:, 2], tseries[:, 3], title=plot_title,
            xlabel="x", ylabel="y", zlabel="z",
            lc=cgrad(:viridis), line_z=time,
            colorbar_title=" \n \ntime",
            right_margin=1.5Plots.mm,
            marker=marker, label=nothing)
        mkpath("Figures/data/")
        savefig(p, "Figures/data/$exp_name.png")
    end

    if SAVE
        dir_path = isempty(save_dir) ? "data/benchmarks/" : save_dir
        mkpath(dir_path)
        npzwrite("data/benchmarks/$exp_name.npy", tseries)
    end
    return tseries
end


"""
    bursting_neuron_regimes(; T=1500, Δt=0.05, transient_T=500, process_noise_level=0.0, PLOT=true, plot_name="", save_dir="", SAVE=true)

creates time series in different bursting neuron regimes for gₙₘ₀ₐ=[3,5,7,9,10,10.2]

Kwargs
------
num_T : number of timesteps, default 15000 (Int)

ΔT : time between steps, default 0.01 (Float)

transient_T : timesteps used as start up to avoid transients, default 2000 (Int)

PLOT : Specifies whether a plot is done, default true (BOOL)

save_dir : directory to save the data in, default data/benchmarks (String)

SAVE : Specifies whether the timeseries is saved or not, default true (BOOL)

process_noise_level : The ratio of process noise, default 0 (Float)

plot_name : the name of the file as which plot and data are saved

Returns
-------
time series : Time series of the specified system (AbstractMatrix)

"""
function bursting_neuron_regimes(; T=1500, Δt=0.05, transient_T=500, process_noise_level=0.0, PLOT=true, plot_name="", save_dir="", SAVE=true)
    μₛ = [3.0, 5.0, 7.0, 9.0, 10.0, 10.2]
    exp_name = isempty(plot_name) ? "bursting_neuron_regimes" : plot_name

    tseries = Vector{AbstractMatrix}()
    for μ in μₛ
        ds = bursting_neuron(u0=[-60.0, 0.0386, 0.0231], gₙₘ₀ₐ=μ)
        ts = generate_trajectories(ds, T, transient_T; Δt, process_noise_level,PLOT=false, model_name="regimes")

        push!(tseries, ts)
        if SAVE
            dir_path = isempty(save_dir) ? "data/benchmarks/" : save_dir
            mkpath(dir_path)
            npzwrite("data/benchmarks/$(exp_name)_$(μ).npy", ts)
        end
    end

    if PLOT
        ps = (plot3d(tseries[i][:, 1], tseries[i][:, 2], tseries[i][:, 3]) for i in 1:6)
        g = plot(ps..., layout=6,
            title=["\n3" "\$g_{nmda}\$ as bifurcation parameter\n5" "\n7" "\n9" "\n10.1" "\n10.2"], titlefont=font(11),
            legend=nothing,
            plot_title="bursting neuron regimes", plot_titlevspan=0.08,
            size=(500, 500))
        mkpath("Figures/data/")
        savefig(g, "Figures/data/$exp_name.png")
    end
end


"""
    ns_3d_benchmark(System::String, ; p_change=linear, T=150, Δt=0.01, transient_T=10, u0=nothing, process_noise_level=0.0, PLOT=true, plot_params=false, plot_name="", plot_title="", snapshots=false, SAVE=true, save_dir="", eval=false, eval_run=0,STD=false)

create a non-stationary benchmark system

Be careful, changing the type of change of parameters might lead to very different trajectories

implemented Systems:

### Lorenz
- ShrinkingLorenz : Starts with the chaotic attractor and shrinks and shifts it a bit; ρ=28->23, σ=10->5, β=8/3->0.5

- PaperLorenzBigChange : The ns system used in Patel et al. 2022 with a quick parameter change

- PaperLorenzSmallChange : The ns system used in Patel et al. 2022 with a slow parameter change

- ExplodingLorenz : Starts with a limit cycle and ends in the 
well known chaotic attractor; ρ=22->28

- ShiftingLorenz : Starts with the chaotic attractor and shifts it "forward"; ρ=28->22

### Bursting Neuron
The timescale is much longer than for the Lorenz, and thus needs a lot more time points

- RampUpBN : Increases bursting behavior, g=>8.2->8.8

StopBurstBN : Starts with the limit cycle, and ends in the standard bursting neuron regime; g=10.1->10.25

Parameters
----------
System : Which system to create (String)

Kwargs
------
p_change : Vector containing the function by which the parameters should change, default [linear, linear, linear]

T :  final time, default 150 (Int)


ΔT : time between steps, default 0.01 (Float)


transient_T : timesteps used as start up to avoid transients, default 2000 (Int)

u0 : Initial condition, default nothing (Vector)


PLOT : Specifies whether a plot is done, default true (BOOL)


save_dir : directory to save the data in, default data/benchmarks (String)


SAVE : Specifies whether the timeseries is saved or not, default true (BOOL)


process_noise_level : The ratio of process noise, default 0 (Float)

plot_name : the name of the file as which plot and data are saved

STD : Return normalized series (Bool)

Returns
-------
time series : Time series of the specified system (AbstractMatrix)
"""
function ns_3d_benchmark(System::String, ; p_change=linear, T=150, Δt=0.01, transient_T=10, u0=nothing, process_noise_level=0.0, PLOT=true, plot_params=false, plot_name="", plot_title="", snapshots=false, SAVE=true, save_dir="", eval=false, eval_run=0,STD=false)
    @assert System in valid_ns_systems "$System is not a valid System: $valid_ns_systems"
    exp_name = isempty(plot_name) ? String(System) : plot_name

    u0 = u0 === nothing ? [0.4, 0.4, 0.8] : u0
    t₀ = 0.0
    if occursin("Paper", System)
        t₀ = -600.0
        T = 800.0
    end

    ns_model, params = ns_benchmark_systems(System, p_change, T; u0, t₀, transient_T)

    pl_params = plot_params ? params : []
    tseries = generate_trajectories(ns_model, T, transient_T; Δt,process_noise_level, model_name=System, PLOT, plot_title, eval, eval_run, pl_params, save_name=exp_name,STD)

    if snapshots
        μ₀ = [params[i](t₀ + transient_T) for i in axes(params, 1)]
        μₑₙ₀ = [params[i](T + transient_T) for i in axes(params, 1)]
        μₛ = [μ₀, μₑₙ₀]
        snap_series = Vector{AbstractMatrix}()
        for (i, μ) in enumerate(μₛ)
            if occursin("Lorenz", System)
                ds = lorenz(p=μ)
            elseif occursin("BN", System)
                ds = bursting_neuron(gₙₘ₀ₐ=μ[1])
            else
                throw("not implemented")
            end
            push!(snap_series, generate_trajectories(ds, T, transient_T; Δt, process_noise_level=0.0, PLOT=false, model_name=System,STD=false))
        end

        comp_plots = plot3d(snap_series[1][:, 1], snap_series[1][:, 2], snap_series[1][:, 3], grid=true, xlabel="x", ylabel="y", zlabel="z", lc=:blue, title="Snapshot comparisons", label="Snapshot at t=0", linealpha=0.8)
        plot3d!(snap_series[2][:, 1], snap_series[2][:, 2], snap_series[2][:, 3], grid=true, lc=:orange, label="Snapshot at t=T")
        mkpath("Figures/snapshots/")
        savefig(comp_plots, "Figures/snapshots/Snapshots_$exp_name.png")

        for (i, series) in enumerate(snap_series)
            std_series = StatsBase.standardize(ZScoreTransform, series, dims=1)
            mkpath("data/snapshots/")
            npzwrite("data/snapshots/snaps_$(exp_name)_$i.npy", std_series)
        end
    end

    if SAVE
        dir_path = isempty(save_dir) ? "data/benchmarks/" : save_dir
        mkpath(dir_path)
        npzwrite("data/benchmarks/$exp_name.npy", tseries)
    end
    return tseries
end


"""
    trial_benchmark(System::String, num_Trials::Int; seq_T=10, Δt=0.01, transient_T=20, plot_title="", PLOT=true, save_dir="", SAVE=true, process_noise_level=0.0, lorenz_sys="")

create a benchmark system in trial form

implemented Systems
- trial_lorenz : a lorenz with parameter shifting from 22->28 across trials
- split_lorenz : a specified ns lorenz split into num_trials

Parameter
---------
System : The name of the benchmark system (String)
num_Trials : Number of trials (Int)

Kwargs
------
num_T : number of timesteps, default 15000 (Int)

ΔT : time between steps, default 0.01 (Float)

transient_T : timesteps used as start up to avoid transients, default 2000 (Int)

plot_title : Title for the plot, default: same as System (String)

PLOT : Specifies whether a plot is done, default true (BOOL)

save_dir : directory to save the data in, default data/benchmarks (String)

SAVE : Specifies whether the timeseries is saved or not (BOOL)

process_noise_level : The ratio of process noise, default 0 (Float)

lorenz_trial_sys : The system used to split the lorenz (String)

Returns
-------
time series : Time series of the specified system (AbstractMatrix)
"""
function trial_benchmark(System::String, num_Trials::Int; seq_T=10, Δt=0.01, transient_T=20, plot_title="", PLOT=true, save_dir="", SAVE=true, process_noise_level=0.0, lorenz_sys="")
    @assert System in valid_trial_systems "$System not in valid_systems:$valid_trial_systems"
    @assert num_Trials > 10 "at least 10 trials need to be in the set"

    if System == "trial_lorenz"
        ts = trial_ns_lorenz(num_Trials, seq_T, Δt, 22.0f0, 28.0f0;process_noise_level, transient_T)
    elseif System == "split_lorenz"
        valid_systems = ["ExplodingLorenz", "ShiftingLorenz", "ShrinkingLorenz", "PaperLorenzBigChange", "PaperLorenzSmallChange"]
        @assert lorenz_sys in valid_systems "$lorenz_sys is not a valid System: $valid_systems"

        ts = split_lorenz_trials(lorenz_sys, num_Trials, seq_T, transient_T, Δt, process_noise_level=process_noise_level)
    end

    if PLOT
        plot_ts = ts
        if num_Trials > 10
            println("only 10 random trials are plotted")
            plot_idx = sort!(sample(2:1:num_Trials-1, 8, replace=false))
            pushfirst!(plot_idx, 1)
            push!(plot_idx, num_Trials)
            plot_ts = plot_ts[plot_idx]
        end


        plots = (plot3d(plot_ts[i][:, 1], plot_ts[i][:, 2], plot_ts[i][:, 3], ticks=false, label=nothing, title="\n\n" * string(plot_idx[i]), titlefont=font(10)) for i in 1:10)
        plot(plots..., layout=grid(2, 5),
            plot_title=plot_title * "\ntrials",
            plot_titlevspan=0.05)
        mkpath("Figures/data/")
        savefig("Figures/data/lorenz_trials.png")
    end

    if SAVE
        dir_path = isempty(save_dir) ? "data/benchmarks/" : save_dir
        mkpath(dir_path)
        num_Trials = System == "trial_lorenz" ? num_Trials += 1 : num_Trials
        ts = reshape(reduce(hcat, ts), num_Trials, Int(seq_T÷Δt)  + 1, 3)
        npzwrite("data/benchmarks/$System.npy", ts)
    end
end