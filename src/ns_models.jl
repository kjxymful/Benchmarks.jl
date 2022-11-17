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
function ns_benchmark_systems(sys::String, par_fun::Function, T::Real; u0=[0.4,0.4,0.8], transient_T=0.0,t₀=0.0)::Tuple{ContinuousDynamicalSystem, AbstractVector}
    if !(sys in valid_ns_systems)
        throw("$sys is not a valid system:$valid_ns_systems")
    end
    tmax = T+t₀
    t_shift = transient_T-t₀
    if occursin("Lorenz", sys)
        σ_init = 10.0
        ρ_init = 35.0
        β_init = 8 / 3
        sys_fun = (du, u, p,t) -> ns_lorenz!(du, u, p, t)
        if sys == "ShrinkingLorenz"
            σ_final = 7.0
            ρ_final = 22.0
            β_final = 1.86
        elseif sys == "ShiftingLorenz"
            σ_final = σ_init
            ρ_final = 22.0
            β_final = β_init
        elseif sys == "ExplodingLorenz"
            σ_final = σ_init
            β_final = β_init
            ρ_final = ρ_init
            ρ_init = 23.0
        elseif sys == "PaperLorenzBigChange"
            σ_final = σ_init
            β_final = β_init
            ρ₀ = 154
            ρ₁ = 8
            τ = 100
            params = [t->linear(σ_init,σ_final,tmax,t-t_shift),t -> exponential(ρ₁,0,0,t-t_shift,offset=ρ₀,τ=τ), t->linear(β_init,β_final,tmax,t-t_shift)]
        elseif sys == "PaperLorenzSmallChange"
            σ_final = σ_init
            β_final = β_init
            ρ₀ = 165.5
            ρ₁ = 0.1
            τ = 100
            params = [t->linear(σ_init,σ_final,tmax,t-t_shift),t -> exponential(ρ₁,0,0,t-t_shift,offset=ρ₀,τ=τ), t->linear(β_init,β_final,tmax,t-t_shift)]
        end
        if !occursin("Paper", sys)
            init_params = [σ_init, ρ_init, β_init]
            final_params = [σ_final, ρ_final, β_final]
        end
    elseif occursin("BN",sys)
        u0 = [-1.33192587,0.00833237,0.01267124]
        sys_fun = (du, u, p,t) -> ns_bursting_neuron!(du, u, p, t)
        if sys == "StopBurstBN"
            g_init = 9.25
            g_final = 10.15
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
    return ContinuousDynamicalSystem(sys_fun, u0, params), params
end