using DynamicalSystems
using Plots
using NPZ: npzwrite
using StatsBase

"""
    std_3d_benchmark(System::String; num_T=15000, ΔT=0.01, transient_T=2000, plot_title="", PLOT=true, save_dir="", SAVE=true, MARKERS=false)

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

Returns
-------
time series : Time series of the specified system (AbstractMatrix)

"""
function std_3d_benchmark(System::String; num_T=15000, ΔT=0.01, transient_T=2000, plot_title="", PLOT=true, save_dir="", SAVE=true, MARKER=false)
    valid_systems = ["standard_bursting", "standard_lorenz", "bursting_limit_cycle", "lorenz_limit_cycle"]
    @assert System in valid_systems "$System is not a valid System: $valid_systems"

    if System == "standard_bursting"
        ds = bursting_neuron()
        μ = 10.2
    elseif System == "standard_lorenz"
        μ = 28.0
        ds = lorenz()
    elseif System == "bursting_limit_cycle"
        ds = bursting_neuron(u0=[-60.0, 0.0386, 0.0231], gₙₘ₀ₐ=10.0)
        μ = 10.0
    elseif System == "lorenz_limit_cycle"
        μ = 24.0
        ds = lorenz()
    end
    tseries = gen_series(ds, num_T, ΔT, transient_T)
    time = 0:ΔT:(num_T*ΔT)

    tseries = StatsBase.standardize(ZScoreTransform, tseries, dims=1)

    if PLOT
        title = isempty(plot_title) ? System : plot_title
        marker = MARKER ? :d : false
        p = plot3d(tseries[:, 1], tseries[:, 2], tseries[:, 3], title=title,
            xlabel="x", ylabel="y", zlabel="z",
            lc=cgrad(:viridis), line_z=time,
            colorbar_title=" \n \ntime",
            right_margin=1.5Plots.mm,
            marker=marker)
        mkpath("Figures/")
        savefig(p, "Figures/$System.png")
    end

    if SAVE
        dir_path = isempty(save_dir) ? "data/benchmarks/" : save_dir
        mkpath(dir_path)
        npzwrite("data/benchmarks/$System.npy", tseries)
    end
    return tseries
end


"""
    bursting_neuron_regimes(; num_T=15000, ΔT=0.01, transient_T=2000, PLOT=true, save_dir="", SAVE=true, MARKER=false)

creates time series in different bursting neuron regimes for gₙₘ₀ₐ=[3,5,7,9,10,10.2]

Kwargs
------
num_T : number of timesteps, default 15000 (Int)

ΔT : time between steps, default 0.01 (Float)

transient_T : timesteps used as start up to avoid transients, default 2000 (Int)

PLOT : Specifies whether a plot is done, default true (BOOL)

save_dir : directory to save the data in, default data/benchmarks (String)

SAVE : Specifies whether the timeseries is saved or not, default true (BOOL)

Returns
-------
time series : Time series of the specified system (AbstractMatrix)

"""
function bursting_neuron_regimes(; num_T=15000, ΔT=0.01, transient_T=2000, PLOT=true, save_dir="", SAVE=true)
    μₛ = [3.0, 5.0, 7.0, 9.0, 10.0, 10.2]

    tseries = Vector{AbstractMatrix}()
    for μ in μₛ
        ds = bursting_neuron(u0=[-60.0, 0.0386, 0.0231], gₙₘ₀ₐ=μ)
        ts = gen_series(ds, num_T, ΔT, transient_T)
        ts = StatsBase.standardize(ZScoreTransform, ts, dims=1)

        push!(tseries, ts)
        if SAVE
            dir_path = isempty(save_dir) ? "data/benchmarks/" : save_dir
            mkpath(dir_path)
            npzwrite("data/benchmarks/bursting_neuron_$(μ).npy", ts)
        end
    end

    if PLOT
        ps = (plot3d(tseries[i][:, 1], tseries[i][:, 2], tseries[i][:, 3]) for i in 1:6)
        g = plot(ps..., layout=6,
            title=["\n3" "\$g_{nmda}\$ as bifurcation parameter\n5" "\n7" "\n9" "\n10" "\n10.2"], titlefont=font(11),
            legend=nothing,
            plot_title="3d bursting neuron regimes", plot_titlevspan=0.08,
            size=(500, 500))
        mkpath("Figures/")
        savefig(g, "Figures/bursting_neuron_regimes.png")
    end
end


"""
    ns_3d_benchmark(System::String, ; p_change=[linear, linear, linear], num_T=15000, ΔT=0.01, transient_T=2000, plot_title="", PLOT=true, SAVE=true)

create a non-stationary benchmark system

Be careful, changing the type of change of parameters might lead to very different trajectories

implemented Systems:
- ExplodingLorenz : Starts with a limit cycle and ends in the 
well known chaotic attractor; ρ=22->28

- ShiftingLorenz : Starts with the chaotic attractor and shifts it "forward"; ρ=28->22

- ShrinkingLorenz : Starts with the chaotic attractor and shrinks and shifts it a bit; ρ=28->23, σ=10->5, β=8/3->0.5


Parameters
----------
System : Which system to create (String)

Kwargs
------
p_change : Vector containing the function by which the parameters should change, default [linear, linear, linear]\n
num_T : number of timesteps, default 15000 (Int)\n

ΔT : time between steps, default 0.01 (Float)\n

transient_T : timesteps used as start up to avoid transients, default 2000 (Int)\n

PLOT : Specifies whether a plot is done, default true (BOOL)\n

save_dir : directory to save the data in, default data/benchmarks (String)\n

SAVE : Specifies whether the timeseries is saved or not, default true (BOOL)\n

Returns
-------
time series : Time series of the specified system (AbstractMatrix)
"""
function ns_3d_benchmark(System::String, ; p_change=[linear, linear, linear], num_T=15000, ΔT=0.01, transient_T=2000, plot_title="", PLOT=true, save_dir="", SAVE=true)
    valid_systems = ["ExplodingLorenz", "ShiftingLorenz", "ShrinkingLorenz"]
    @assert System in valid_systems "$System is not a valid System: $valid_systems"

    tmax = (num_T + transient_T) * ΔT

    ns_model = ns_lorenz_systems(System, p_change, tmax)

    tseries = generate_ns_trajectories(ns_model, tmax, transient_T, Δt=ΔT, PLOT=PLOT)

    tseries = StatsBase.standardize(ZScoreTransform, tseries, dims=1)

    if SAVE
        dir_path = isempty(save_dir) ? "data/benchmarks/" : save_dir
        mkpath(dir_path)
        npzwrite("data/benchmarks/$System.npy", tseries)
    end
    return tseries
end

