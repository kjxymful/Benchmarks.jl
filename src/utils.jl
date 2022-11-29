using Plots
using DifferentialEquations: Tsit5, DiscreteCallback

# adjust the correct tipping points
TP_dict = Dict("StopBurstBN" => 10.0, 
                "RampUpBN"=>2.6, 
                "ShrinkingLorenz"=>24, 
                "ShiftingLorenz"=>24, 
                "ExplodingLorenz"=>24,
                "PaperLorenzBigChange"=>166,
                "PaperLorenzSmallChange"=>166)
TP_loc = Dict("StopBurstBN" => 25001, 
                "RampUpBN"=>9001, 
                "ShrinkingLorenz"=>12694, 
                "ShiftingLorenz"=>12694, 
                "ExplodingLorenz"=>1251,
                "PaperLorenzBigChange"=>64056,
                "PaperLorenzSmallChange"=>76096)


function plot_param(params, time, Ttr,model_name,TP)
    time1 = time .+ Ttr
    if occursin("Lorenz",model_name)
        p = plot(time, params[1].(time1)/params[1](time1[1]), label="σ₀=$(round(params[1](time1[1]),digits=1))", legend=:inside, ylabel="% of IC")
        plot!(time, params[2].(time1)/params[2](time1[1]), label="ρ₀=$(round(params[2](time1[1]),digits=1))")
        plot!(time, params[3].(time1)/params[3](time1[1]), label="β₀=$(round(params[3](time1[1]),digits=1))")
        μ = params[2]
    elseif occursin("BN",model_name)
        p = plot(time, params[1].(time1)/params[1](time1[1]), label="gₙₘ₀ₐ₀=$(round(params[1](time1[1]),digits=1))",legend=:inside, ylabel="% of IC")
        μ = params[1]
    else
        throw("No parameter specifices for this model implemented")
    end
    if TP
        TP_val = TP_dict[model_name]
        # add a vertical line at time where params[1] equals TP_val
        if μ(0) > μ(10)
            TP_location = findfirst(x->x<=TP_val, μ.(time1))
        elseif μ(0) < μ(10)
            TP_location = findfirst(x->x>=TP_val, μ.(time1))
        else
            throw("Something is weird") 
        end           
        @show μ.(time1)[1], μ.(time1)[end], TP_location
        vline!(p,[time[TP_location]], label="TP at :$(Int(round(TP_location,digits=1)))")
    end    
    return p
end
    

"""
    generate_trajectories(model, T, transient_T; Δt=0.01, process_noise_level=0.0, STD=true, save_name="", PLOT=true, plot_title="", pl_params=zeros(3), eval=false, eval_run=0, model_name="")::AbstractMatrix

generate the trajectory of a given ns system

Parameters
----------
ns_model : an initalized dynamical System
t₀ : Start time, default 0 (Float64)
STD : true for a standardized time series
PLOT : If a plot should be done, default true (BOOL)
plot_title : Specifies plot title
save_name : as what to save the file

*Plotting to Figures/*
"""
function generate_trajectories(model::GeneralizedDynamicalSystem, T::Real, transient_T::Real; Δt=0.01, process_noise_level=0.0, STD=true, save_name="", PLOT=true, plot_title="", pl_params=zeros(3), eval=false, eval_run=0, model_name="",TP=false)::AbstractMatrix
    std_ = process_noise_level==0 ? [0,0,0] : std_model(model, T, transient_T;Δt)
    ts = trajectory(model, T; Δt, Ttr=transient_T, diffeq=(alg=Tsit5(), callback=dynamical_noise_callback(process_noise_level, std_)))
    u = Matrix(ts)
    u_std = StatsBase.standardize(ZScoreTransform, u, dims=1)

    t = 0.0:Δt:T
    if PLOT
        if occursin("Paper", model_name)
            uS0 = u_std[findall(x -> x <= 600, t), :]
            uB0 = u_std[findall(x -> x > 600, t), :]
            p = plot3d(uS0[:, 1], uS0[:, 2], uS0[:, 3],
                grid=true,
                lc=:red, label="t\$\\leq\$0 (600)", linealpha=1, title=plot_title)
            plot3d!(p, uB0[:, 1], uB0[:, 2], uB0[:, 3], lc=:black, label="t>0 (600)", linealpha=0.8)
        else
            p = plot3d(u_std[1:end, 1], u_std[1:end, 2], u_std[1:end, 3],
                grid=true,
                lc=cgrad(:viridis), line_z=t[1:end],
                colorbar_title=" \n \ntime",
                right_margin=1.5Plots.mm, title=plot_title,label=nothing)
        end
        if pl_params != zeros(3)
            p_par = plot_param(pl_params, t, transient_T, model_name,TP)
            l = grid(2,1, heights=[0.8,0.2])
            p = plot(p, p_par, layout=l, plot_title=plot_title)
        end
        if eval
            mkpath("Figures/eval_$(model_name)")
            savefig(p, "Figures/eval_$(model_name)/$save_name $eval_run.png")
        else
            save_name = isempty(save_name) ? model_name : save_name
            mkpath("Figures/data/")
            savefig(p, "Figures/data/$save_name.png")
        end
    end
    if STD 
        return u_std
    else
        return u
    end
end



function sigmoid(lower::Real, upper::Real, tmax::Real, x::Real)
    return lower + (upper - lower) / (1 + exp(-2.5 * x + tmax))
end

function linear(start::Real, end_::Real, tmax::Real, x::Real)
    return start + (end_ - start) / tmax * x
end

function exponential(start::Real, end_::Real, tmax::Real, x::Real; offset=0.0, τ=0.0)
    if τ == 0
        τ = tmax / log(end_ / start)
    end
    return offset + start * exp(x / τ)
end


function check_validity(sys::String)
    if !any(x->sys in x, [valid_ns_systems,valid_trial_systems, valid_std_systems,valid_regimes])
        throw("$benchmark_system is not a valid system")
    end
    return true
end

function dynamical_noise_callback(dynamical_noise_level::AbstractFloat, std::AbstractVector)::DiscreteCallback
    noise_level = dynamical_noise_level .* std
    dn = dynamical_noise_level == 0 ? false : true
    condition(u, t, integrator) = dn
    affect!(integrator) = set_state!(integrator, get_state(integrator) .+ ϵ.(noise_level))
    cb = DiscreteCallback(condition, affect!)
    return cb
end

function std_model(ds::ContinuousDynamicalSystem,T::Real, transient_T::Real; Δt=0.01)
    tseries = Matrix(trajectory(ds, T, Ttr=transient_T,Δt=Δt))
    std_ = [std(tseries[:, i]) for i in axes(tseries, 2)]
    return std_
end

# function std_model(ds::ODEProblem, tmax,transient_T;Δt=0.01,t₀=0.0)
#     cde = ContinuousDynamicalSystem(ds)
#     return std_model(cde, tmax,transient_T;Δt,t₀)
# end

ϵ(process_noise::AbstractFloat)::AbstractFloat = randn() * process_noise
