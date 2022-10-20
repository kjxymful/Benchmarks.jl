using StatsBase

"""
an object defining the pre chosen ns systems

params need to be funtions f(t)
"""
mutable struct ns_systems{V<:AbstractVector,Vf<:AbstractVector}
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
function ns_lorenz_systems(sys::String, par_fun::AbstractVector, tmax::AbstractFloat; u0=[0.5, 0.5, 0.5])
    valid_models = ["ExplodingLorenz", "ShrinkingLorenz", "ShiftingLorenz"]
    if !(sys in valid_models)
        throw("$sys is not a valid system")
    end
    if occursin("Lorenz", sys)
        σ_init = 10.0
        ρ_init = 28.0
        β_init = 8 / 3
        sys_fun = ns_lorenz!
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
        end
        init_params = [σ_init, ρ_init, β_init]
        final_params = [σ_final, ρ_final, β_final]
    else
        init_params = zeros(3)
        final_params = zeros(3)
        throw("not implemented")
    end
    params = [t -> par_fun[i](init_params[i], final_params[i], tmax, t) for i in axes(init_params, 1)]
    return ns_systems(sys_fun, params, u0, sys)
end


