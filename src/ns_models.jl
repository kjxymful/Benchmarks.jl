using StatsBase

"""
an object defining the pre chosen ns systems

params need to be funtions f(t)
"""
struct ns_systems{V<:AbstractVector,Vf<:AbstractVector} <: AbstractDynamicalSystem
    sys::Function
    params::Vf
    u0::V
    name::String
end

"""
choose and initialize the given ns system

Parameters
----------
sys : The type of Benchmark system (String)\n
par_fun : the function according to which the parameters will evolve (Vector{Function})

Returns
-------
Initialized ns_systems for the lorenz ns benchmarks
"""
function ns_benchmark_systems(sys::String, par_fun::Function, tmax::AbstractFloat, process_noise_level::Vector{Float64}; u0=[0.5, 0.5, 0.5], transient_T=0.0,t₀=0.0)
    valid_models = ["RampUpBN", "SuddenBurstBN", "ExplodingLorenz", "ShrinkingLorenz", "ShiftingLorenz","PaperLorenzBigChange", "PaperLorenzSmallChange"]
    if !(sys in valid_models)
        throw("$sys is not a valid system:$valid_models")
    end
    if occursin("Lorenz", sys)
        σ_init = 10.0
        ρ_init = 28.0
        β_init = 8 / 3
        sys_fun = (du, u, p,t) -> ns_lorenz!(du, u, p, t, process_noise=process_noise_level)
        if sys == "ShrinkingLorenz"
            σ_final = 5.0
            ρ_final = 23.0
            β_final = 0.5
        elseif sys == "ShiftingLorenz"
            σ_final = σ_init
            ρ_final = 22.0
            β_final = β_init
        elseif sys == "ExplodingLorenz"
            σ_final = σ_init
            β_final = β_init
            ρ_final = ρ_init
            ρ_init = 22.0
        elseif sys == "PaperLorenzBigChange"
            σ_final = σ_init
            β_final = β_init
            ρ₀ = 154
            ρ₁ = 8
            τ = 100
            params = [t->linear(σ_init,σ_final,tmax+transient_T,t-transient_T),t -> exponential(ρ₁,0,0,t-transient_T,offset=ρ₀,τ=τ), t->linear(β_init,β_final,tmax+transient_T,t-transient_T)]
        elseif sys == "PaperLorenzSmallChange"
            σ_final = σ_init
            β_final = β_init
            ρ₀ = 165.5
            ρ₁ = 0.1
            τ = 100
            params = [t->linear(σ_init,σ_final,tmax+transient_T,t-transient_T),t -> exponential(ρ₁,0,0,t-transient_T,offset=ρ₀,τ=τ), t->linear(β_init,β_final,tmax+transient_T,t-transient_T)]
        end
        if !occursin("Paper", sys)
            init_params = [σ_init, ρ_init, β_init]
            final_params = [σ_final, ρ_final, β_final]
        end
    elseif occursin("BN",sys)
        u0 = [-20,0.5,0.03]
        sys_fun = (du, u, p,t) -> ns_bursting_neuron!(du, u, p, t, process_noise=process_noise_level)
        if sys == "SuddenBurstBN"
            g_init = 10.08
            g_final = 10.4
        elseif sys == "RampUpBN"
            g_init = 2
            g_final = 4
        end
        init_params = [g_init]
        final_params = [g_final]  
    end
    if !occursin("Paper",sys)
        params = [t -> par_fun(init_params[i], final_params[i], tmax, t-transient_T) for i in axes(init_params, 1)]
    end
    return ContinuousDynamicalSystem(sys_fun, u0, params,t0=t₀), params
end


function ns_systems_bench(sys::String, par_fun::Function, tmax::AbstractFloat, process_noise_level::Vector{Float64}; u0=[0.5, 0.5, 0.5], transient_T=0.0,t₀=0.0)
    valid_models = ["RampUpBN", "SuddenBurstBN", "ExplodingLorenz", "ShrinkingLorenz", "ShiftingLorenz","PaperLorenzBigChange", "PaperLorenzSmallChange"]
    if !(sys in valid_models)
        throw("$sys is not a valid system:$valid_models")
    end
    if occursin("Lorenz", sys)
        σ_init = 10.0
        ρ_init = 28.0
        β_init = 8 / 3
        sys_fun = (du, u, p,t) -> ns_lorenz!(du, u, p, t, process_noise=process_noise_level)
        if sys == "ShrinkingLorenz"
            σ_final = 5.0
            ρ_final = 23.0
            β_final = 0.5
        elseif sys == "ShiftingLorenz"
            σ_final = σ_init
            ρ_final = 22.0
            β_final = β_init
        elseif sys == "ExplodingLorenz"
            σ_final = σ_init
            β_final = β_init
            ρ_final = ρ_init
            ρ_init = 22.0
        elseif sys == "PaperLorenzBigChange"
            σ_final = σ_init
            β_final = β_init
            ρ₀ = 154
            ρ₁ = 8
            τ = 100
            params = [t->linear(σ_init,σ_final,tmax,t-transient_T),t -> exponential(ρ₁,0,0,t-transient_T,offset=ρ₀,τ=τ), t->linear(β_init,β_final,tmax,t-transient_T)]
        elseif sys == "PaperLorenzSmallChange"
            σ_final = σ_init
            β_final = β_init
            ρ₀ = 165.5
            ρ₁ = 0.1
            τ = 100
            params = [t->linear(σ_init,σ_final,tmax,t-transient_T),t -> exponential(ρ₁,0,0,t-transient_T,offset=ρ₀,τ=τ), t->linear(β_init,β_final,tmax,t-transient_T)]
        end
        if !occursin("Paper", sys)
            init_params = [σ_init, ρ_init, β_init]
            final_params = [σ_final, ρ_final, β_final]
        end
    elseif occursin("BN",sys)
        u0 = [-20,0.5,0.03]
        sys_fun = (du, u, p,t) -> ns_bursting_neuron!(du, u, p, t, process_noise=process_noise_level)
        if sys == "SuddenBurstBN"
            g_init = 10.08
            g_final = 10.4
        elseif sys == "RampUpBN"
            g_init = 2
            g_final = 4
        end
        init_params = [g_init]
        final_params = [g_final]  
    end
    if !occursin("Paper",sys)
        params = [t -> par_fun(init_params[i], final_params[i], tmax, t-transient_T) for i in axes(init_params, 1)]
    end
    return ns_systems(sys_fun, params,u0,sys)
end

