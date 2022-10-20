using DifferentialEquations

function gen_series(ds::GeneralizedDynamicalSystem, num_T::Int, ΔT::AbstractFloat, transient_T::Int)
    T_end = num_T * ΔT
    ts = trajectory(ds, T_end; Δt=ΔT, Ttr=transient_T)
    return Matrix(ts)
end


"""
    generate_ns_trajectories(ns_model::ns_systems,tmax::AbstractFloat, skip_steps::Int;Δt=0.01, t₀=0.0, PLOT=true)

generate the trajectory of a given ns system

Parameters
----------
ns_model : an initalized ns system (ns_systems)
t₀ : Start time, default 0 (Float64)
PLOT : If a plot should be done, default true (BOOL)

*Plotting to Figures/*
"""
function generate_ns_trajectories(ns_model::ns_systems,tmax::AbstractFloat, skip_steps::Int;Δt=0.01, t₀=0.0, PLOT=true)
    tspan = (t₀, tmax)
    prob = ODEProblem(ns_model.sys, ns_model.u0, tspan, ns_model.params)
    sol = solve(prob)

    u = [sol(t) for t in (t₀+skip_steps*Δt):Δt:tmax]
    u = permutedims(hcat(u...))
    u = StatsBase.standardize(ZScoreTransform, u, dims=1)
    t = collect(t₀:0.01:(tmax-skip_steps*Δt))

    if PLOT
        p = plot3d(u[1:end, 1], u[1:end, 2], u[1:end, 3],
            xlabel="x", ylabel="y", zlabel="z",
            lc=cgrad(:viridis), line_z=t[1:end] * 100,
            colorbar_title=" \n \ntime",
            right_margin=1.5Plots.mm)
        mkpath("Figures/")
        savefig(p, "Figures/$(ns_model.name).png")
    end
    return u
end


function sigmoid(lower::Real, upper::Real, tmax::Real, x::Real)
    return lower + (upper - lower) / (1 + exp(-2.5 * x + tmax))
end

function linear(start::Real, end_::Real, tmax::Real, x::Real)
    return start + (end_ - start) / tmax * x
end

function exponential(start::Real, end_::Real, tmax::Real, x::Real)
    τ = tmax/log(end_/start)
    return start * exp(x / τ)
end