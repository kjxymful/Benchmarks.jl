
abstract type AbstractDynamicalSystem end

struct ds_sytem{V<:AbstractVector} <: AbstractDynamicalSystem
    sys::Function
    params::V
    u0::V
    name::String
end



"""
creates the bursting neuron model by Durstewitz 2009
as ContinuousDynamicalSystem

Parameters:
-----------
gₙₘ₀ₐ (Float): the bifurcation parameter
u0 (Vector{Float}): initial condition

every other parameter specified accessible as in Paper\n 
(with subscripts)

Returns:
--------
ds (ContinuousDynamicalSystem): dynamical system with the bursting neuron ODES
...
"""
function bursting_neuron(; u0=[-24.4694, 0.0386, 0.0231],
    I=0,
    Cₘ=6,
    gₗ=8,
    Eₗ=-80,
    gₙₐ=20,
    Eₙₐ=60,
    Vₕₙₐ=-20,
    kₙₐ=15,
    gₖ=10,
    Eₖ=-90,
    Vₕₖ=-25,
    kₖ=5,
    τₙ=1,
    gₘ=25,
    Vₕₘ=-15,
    kₘ=5,
    τₕ=200,
    gₙₘ₀ₐ=10.2,t₀=0.0)::ContinuousDynamicalSystem
    p = [I, Cₘ, gₗ, Eₗ, gₙₐ, Eₙₐ, Vₕₙₐ, kₙₐ, gₖ, Eₖ, Vₕₖ, kₖ, τₙ, gₘ, Vₕₘ, kₘ, τₕ, gₙₘ₀ₐ]
    rhs = (du,u,par,t) -> loop_burstn!(du,u,par,t)
    ds = ContinuousDynamicalSystem(rhs, u0, p,t0=t₀)
    return ds
end

"""
uses DifferentialEquations somehow ends up with different attractors
"""
function bursting_neuron(method_holder::String; u0=[-24.4694, 0.0386, 0.0231],
    I=0.0,
    Cₘ=6.0,
    gₗ=8.0,
    Eₗ=-80.0,
    gₙₐ=20.0,
    Eₙₐ=60.0,
    Vₕₙₐ=-20.0,
    kₙₐ=15.0,
    gₖ=10.0,
    Eₖ=-90.0,
    Vₕₖ=-25.0,
    kₖ=5.0,
    τₙ=1.0,
    gₘ=25.0,
    Vₕₘ=-15.0,
    kₘ=5.0,
    τₕ=200.0,
    gₙₘ₀ₐ=10.2)::ds_sytem
    p = [I, Cₘ, gₗ, Eₗ, gₙₐ, Eₙₐ, Vₕₙₐ, kₙₐ, gₖ, Eₖ, Vₕₖ, kₖ, τₙ, gₘ, Vₕₘ, kₘ, τₕ, gₙₘ₀ₐ]
    ds = ds_sytem((du, u, param, t) -> loop_burstn!(du, u, param, t), p, u0, "bursting_neuron")
    return ds
end


function lorenz(;u0=[0.5,0.5,0.5], ρ=28.0, σ=10.0, β=8/3, p=[],t₀=0.0)::ContinuousDynamicalSystem
    if isempty(p)
        p = [σ, ρ, β]
    else
        p=p
    end
    rhs = (du,u,par,t) -> loop_lorenz!(du,u,par,t)
    return ContinuousDynamicalSystem(rhs, u0,p, Systems.lorenz_iip().jacobian,t0=t₀)
end


function lorenz(method_holder::String; u0=[0.5,0.5,0.5], ρ=28.0, σ=10.0, β=8/3, p=[])::ds_sytem
    if isempty(p)
        p = [σ, ρ, β]
    else
        p=p
    end
    return ds_sytem((du,u,p,t)->loop_lorenz!(du,u,p,t), p, u0, "lorenz")
end