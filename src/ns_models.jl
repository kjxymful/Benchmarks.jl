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
function ns_lorenz_systems(sys::String, par_fun::AbstractVector, tmax::AbstractFloat, process_noise_level::Vector{Float64}; u0=[0.5, 0.5, 0.5])
    valid_models = ["ExplodingLorenz", "ShrinkingLorenz", "ShiftingLorenz","PaperLorenzBigChange", "PaperLorenzSmallChange"]
    if !(sys in valid_models)
        throw("$sys is not a valid system:$valid_models")
    end
    if occursin("Lorenz", sys)
        σ_init = 10.0
        ρ_init = 28.0
        β_init = 8 / 3
        sys_fun(du, u, p,t) = ns_lorenz!(du, u, p, t, process_noise=process_noise_level)
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
            params = [t->linear(σ_init,σ_final,tmax,t),t -> exponential(ρ₁,0,0,t,offset=ρ₀,τ=τ), t->linear(β_init,β_final,tmax,t)]
        elseif sys == "PaperLorenzSmallChange"
            σ_final = σ_init
            β_final = β_init
            ρ₀ = 165.5
            ρ₁ = 0.1
            τ = 100
            params = [t->linear(σ_init,σ_final,tmax,t),t -> exponential(ρ₁,0,0,t,offset=ρ₀,τ=τ), t->linear(β_init,β_final,tmax,t)]
        end
        if !occursin("Paper", sys)
            init_params = [σ_init, ρ_init, β_init]
            final_params = [σ_final, ρ_final, β_final]
        end
    else
        init_params = zeros(3)
        final_params = zeros(3)
        throw("not implemented")
    end
    if !occursin("Paper",sys)
        params = [t -> par_fun[i](init_params[i], final_params[i], tmax, t) for i in axes(init_params, 1)]
    end
    return ns_systems(sys_fun, params, u0, sys)
end


