using DifferentialEquations

function generate_trajetories(ds::GeneralizedDynamicalSystem, num_T::Int, ΔT::AbstractFloat, transient_T::Int)
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
function generate_trajectories(model::AbstractDynamicalSystem,tmax::AbstractFloat, skip_steps::Int;Δt=0.01, t₀=0.0, PLOT=true, plot_title="", eval=false, eval_run=0, save_name="")
    tspan = (t₀, tmax)

    prob = ODEProblem(model.sys, model.u0, tspan, model.params)
    sol = solve(prob)

    u = [sol(t) for t in (t₀+skip_steps*Δt):Δt:tmax]
    u = permutedims(hcat(u...))
    t = collect(t₀:0.01:(tmax-skip_steps*Δt))
    u_std = StatsBase.standardize(ZScoreTransform, u, dims=1)

    if PLOT
        if occursin("Paper", model.name)
            uS0 = u_std[findall(x->x<=0,t),:]
            uB0 = u_std[findall(x->x>0,t),:]
            p = plot3d(uS0[:, 1], uS0[:, 2], uS0[:, 3],
                xlabel="x", ylabel="y", zlabel="z",
                lc=:red,label="t\$\\leq\$0",linealpha=1, title=plot_title)
            plot3d!(p, uB0[:,1],uB0[:,2],uB0[:,3],lc=:black,label="t>0",linealpha=0.8)
        else
            p = plot3d(u_std[1:end, 1], u_std[1:end, 2], u_std[1:end, 3],
                xlabel="x", ylabel="y", zlabel="z",
                lc=cgrad(:viridis), line_z=t[1:end] * 100,
                colorbar_title=" \n \ntime",
                right_margin=1.5Plots.mm, title=plot_title)
        end
        if eval
            mkpath("Figures/eval_$(model.name)")
            savefig(p, "Figures/eval_$(model.name)/$eval_run.png")
        else
            save_path = save_name == "" ? "Figures/$(model.name).png" : "Figures/"*save_name*".png"
            mkpath("Figures/")
            savefig(p, save_path)
        end
    end
    return u
end


function sigmoid(lower::Real, upper::Real, tmax::Real, x::Real)
    return lower + (upper - lower) / (1 + exp(-2.5 * x + tmax))
end

function linear(start::Real, end_::Real, tmax::Real, x::Real)
    return start + (end_ - start) / tmax * x
end

function exponential(start::Real, end_::Real, tmax::Real, x::Real;offset=0.0,τ=0.0)
    if τ == 0
        τ = tmax/log(end_/start)
    end
    return offset+start * exp(x / τ)
end

