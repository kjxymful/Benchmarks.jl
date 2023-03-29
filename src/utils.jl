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
model: An instance of ContinuousDynamicalSystem
T: The time to simulate the system for (Real)
transient_T: The time to simulate the system for before saving the trajectory (Real)
Δt: The time step to use for the simulation (Real)
process_noise_level: The level of process noise to add to the system (Real)
STD: Whether to standardize the trajectory (Bool)
save_name: The name of the file to save the trajectory to (String)
PLOT: Whether to plot the trajectory (Bool)
plot_title: The title of the plot (String)


*Plotting to Figures/*
"""
function generate_trajectories(model::ContinuousDynamicalSystem, T::Real, transient_T::Real; Δt=0.01, process_noise_level=0.0, STD=true, save_name="", PLOT=true, plot_title="", pl_params=zeros(3), eval=false, eval_run=0, model_name="",TP=false)::AbstractMatrix
    std_ = process_noise_level==0 ? [0,0,0] : std_model(model, T, transient_T;Δt)
    ts = trajectory(model, T; Δt, Ttr=transient_T, diffeq=(alg=Tsit5(), callback=dynamical_noise_callback(process_noise_level, std_)))
    u = Matrix(ts)
    u_std = StatsBase.standardize(ZScoreTransform, u, dims=1)

    t = 0.0:Δt:T
    if PLOT
        if occursin("Paper", model_name) && !occursin("complicated", save_name)
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

function complicated(start::Real, end_::Real, tmax::Real, x::Real)
    e_val = end_ - start
    return (abs(x / tmax * cos(x * 2 * pi / tmax))) * e_val + start
end

function complicated_lorenz(start::Real, end_::Real, tmax::Real, x::Real; offset=0.0, τ=0.0)
    s_val = offset == 0 ? start : offset
    e_val = τ == 0 ? end_ : exponential(start, end_, tmax, 200; offset=offset, τ=τ)
    bump(x) = (x+600)/800
    e_val = e_val - bump(800) - s_val
    return (abs((x+600)/800*sin(x*2*pi/800)))*e_val+bump(x)+ s_val
end


function check_validity(sys::String)
    if !any(x->sys in x, [valid_ns_systems,valid_trial_systems, valid_std_systems,valid_regimes])
        throw("$sys is not a valid system")
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

ϵ(process_noise::AbstractFloat)::AbstractFloat = randn() * process_noise

"""
Snapshot generator for Paper Lorenz systems
"""
function lorenz_at_t(T::AbstractFloat;p_fun="", STD=true)
    ρ₀ = 154
    ρ₁ = 8
    τ = 100
    t₀ = -600
    Δt = 0.01
    tmax = 200
    time = -600:0.01:200
    T = time[Int(T*80000)+1]
    p_fun = isempty(p_fun) ? "exponential" : p_fun
    p_sym = Symbol(p_fun)
    p_change = @eval $p_sym
    ρ(t) = p_change(ρ₁,0,0,t,offset=ρ₀,τ=τ)
    ρ = ρ(T)
    ds = lorenz(;ρ)
    ts = Matrix(trajectory(ds, 800,Ttr=20))
    if STD
        ns_sys,_ = ns_benchmark_systems("PaperLorenzBigChange", p_change, tmax; t₀, transient_T=20)
        nsts = generate_trajectories(ns_sys, 800, 20, PLOT=false, STD=false)
        std_ = std(nsts, dims=1)
        mean_ = mean(nsts,dims=1)
        for i in 1:3
            global ts[:,i] = standardize_to_nsseries(ts[:,i], mean_[i], std_[i])
        end
    end
    return ts
end
    

"""
Snapshot generator for Bursting Neuron Systems
"""
function BN_at_t(T::AbstractFloat;p_fun="",STD=true,pn=false,sys="StopBurstBN")
    g_init = 9.25
    g_final = 10.15
    if sys == "RampUpBN"
        g_init = 8.2
        g_final = 8.8
    elseif sys == "EasyBN"
        g_init = 3
        g_final = 4.5
    elseif sys == "SimpleRampBN"
        g_init = 7.0
        g_final = 8.8
    end
    tmax = 1500
    Δt = 0.05
    time = 0:Δt:tmax
    t_id = Int(T*tmax/0.05)+1
    T = time[t_id]
    p_fun = isempty(p_fun) ? "linear" : p_fun
    p_sym = Symbol(p_fun)
    p_change = @eval $p_sym
    ns_sys,_ = ns_benchmark_systems(sys, p_change, tmax; transient_T=500)
    nsts = generate_trajectories(ns_sys, 1500, 20, PLOT=false, STD=false, Δt=Δt)
    g(t) = p_change(g_init, g_final, tmax, t)
    gₙₘ₀ₐ = g(T)
    ds = bursting_neuron(;gₙₘ₀ₐ)
    kwargs = ()
    if pn
        kwargs = (diffeq=(alg=Tsit5(), callback=dynamical_noise_callback(0.05, std(nsts,dims=1))))
    end
    ts = Matrix(trajectory(ds, tmax; Ttr=500, Δt, kwargs...))
    if STD
        std_ = std(nsts,dims=1)
        mean_ = mean(nsts,dims=1)
        for i in 1:3
            ts[:,i] = standardize_to_nsseries(ts[:,i], mean_[i], std_[i])
        end
    end
    return ts
end

"""
Snapshot generator for Roessler
"""
function rssl_at_t(T::AbstractFloat; p_fun="", STD=true)
    init_params = [0.1, 0.1, 14.0f0]
    final_params = [0.2, 0.2, 5.7]
    tmax = 150
    Δt = 0.01
    TTr = 20
    time = 0f0:Δt:tmax
    t_id = Int(T * tmax / 0.01) + 1
    T = time[t_id]
    p_fun = isempty(p_fun) ? "linear" : p_fun
    p_sym = Symbol(p_fun)
    p_change = @eval $p_sym
    p(t) = p_change.(init_params, final_params, Ref(tmax), Ref(t))
    p = p(T)
    a,b,c= p
    ds = Systems.roessler(;a,b,c)
    ts = Matrix(trajectory(ds, tmax, Ttr=20))
    if STD
        ns_sys, _ = ns_benchmark_systems("Roessler", p_change, tmax; transient_T=20)
        nsts = generate_trajectories(ns_sys, tmax, TTr, PLOT=false, STD=false)
        std_ = std(nsts, dims=1)
        mean_ = mean(nsts, dims=1)
        for i in 1:3
            global ts[:, i] = standardize_to_nsseries(ts[:, i], mean_[i], std_[i])
        end
    end
    return ts
end

function standardize_to_nsseries(u::AbstractVector, mean::AbstractFloat, std::AbstractFloat)
    return (u .- mean) ./ std
end