using StatsBase
using DynamicalSystems

abstract type AbstractDynamicalSystem end

struct ds_sytem{V<:AbstractVector} <: AbstractDynamicalSystem
    sys::Function
    params::V
    u0::V
    name::String
end


@inline @inbounds function loop_burstn!(du,u, p, t;process_noise=zeros(3))
    I = p[1]
    Cₘ = p[2]
    gₗ = p[3]
    Eₗ = p[4]
    gₙₐ = p[5]
    Eₙₐ = p[6]
    Vₕₙₐ = p[7]
    kₙₐ = p[8]
    gₖ = p[9]
    Eₖ = p[10]
    Vₕₖ = p[11]
    kₖ = p[12]
    τₙ = p[13]
    gₘ = p[14]
    Vₕₘ = p[15]
    kₘ = p[16]
    τₕ = p[17]
    gₙₘ₀ₐ = p[18]
    Eₙₘ₀ₐ = 0 # as far as i could find

    V = u[1]
    n = u[2]
    h = u[3]
    s∞(V) = 1 / (1 + 0.33 * exp(-0.0625 * V))
    m∞(V, Vₕₙₐ, kₙₐ) = 1 / (1 + exp((Vₕₙₐ - V) / kₙₐ))
    n∞(V, Vₕₖ, kₖ) = 1 / (1 + exp((Vₕₖ - V) / kₖ))
    h∞(V, Vₕₘ, kₘ) = 1 / (1 + exp((Vₕₘ - V) / kₘ))

    du[1] = (I - gₗ * (V - Eₗ) - gₙₐ * m∞(V, Vₕₙₐ, kₙₐ) * (V - Eₙₐ)
           -
           gₖ * n * (V - Eₖ) - gₘ * h * (V - Eₖ)
           -
           gₙₘ₀ₐ * s∞(V) * (V - Eₙₘ₀ₐ)) / Cₘ + ϵ(process_noise[1])
    du[2] = (n∞(V, Vₕₖ, kₖ) - n) / τₙ + ϵ(process_noise[2])
    du[3] = (h∞(V, Vₕₘ, kₘ) - h) / τₕ + ϵ(process_noise[3])
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
    gₙₘ₀ₐ=10.2,
    process_noise=zeros(3))
    p = [I, Cₘ, gₗ, Eₗ, gₙₐ, Eₙₐ, Vₕₙₐ, kₙₐ, gₖ, Eₖ, Vₕₖ, kₖ, τₙ, gₘ, Vₕₘ, kₘ, τₕ, gₙₘ₀ₐ]
    ds = ds_sytem((du,u,p,t)->loop_burstn!(du,u,p,t,process_noise=process_noise), p, u0, "bursting_neuron")
    return ds
end

@inline @inbounds function loop_lorenz!(du,u,p,t; process_noise=zeros(3))
    σ, ρ, β = p
    du[1] = σ * (u[2] - u[1]) + ϵ(process_noise[1])
    du[2] = u[1] * (ρ - u[3]) - u[2] + ϵ(process_noise[2])
    du[3] = u[1] * u[2] - β * u[3] + ϵ(process_noise[3])
end

function lorenz(; u0=[0.5,0.5,0.5], ρ=28.0, σ=10.0, β=8/3, process_noise=zeros(3))
    p = [σ, ρ, β]
    return ds_sytem((du,u,p,t)->loop_lorenz!(du,u,p,t,process_noise=process_noise), p, u0, "lorenz")
end


function ns_lorenz!(du, u, p, t; process_noise=zeros(3))
    σ, ρ, β = p
    du[1] = σ(t) * (u[2] - u[1]) + ϵ(process_noise[1])
    du[2] = u[1] * (ρ(t) - u[3]) - u[2] + ϵ(process_noise[2])
    du[3] = u[1] * u[2] - β(t) * u[3] + ϵ(process_noise[3])
end

ϵ(process_noise::AbstractFloat) = randn()*process_noise
