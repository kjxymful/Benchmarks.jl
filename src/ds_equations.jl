using StatsBase


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
    Eₙₘ₀ₐ = 0

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

@inline @inbounds function loop_lorenz!(du, u, p, t; process_noise=zeros(3))
    σ, ρ, β = p
    du[1] = σ * (u[2] - u[1]) + ϵ(process_noise[1])
    du[2] = u[1] * (ρ - u[3]) - u[2] + ϵ(process_noise[2])
    du[3] = u[1] * u[2] - β * u[3] + ϵ(process_noise[3])
end

@inline @inbounds function ns_lorenz!(du, u, p, t; process_noise=zeros(3))
    σ, ρ, β = p
    du[1] = σ(t) * (u[2] - u[1]) + ϵ(process_noise[1])
    du[2] = u[1] * (ρ(t) - u[3]) - u[2] + ϵ(process_noise[2])
    du[3] = u[1] * u[2] - β(t) * u[3] + ϵ(process_noise[3])
end

@inline @inbounds function ns_bursting_neuron!(du,u, p, t;process_noise=zeros(3))
    I = 0
    Cₘ = 6
    gₗ = 8
    Eₗ = -80
    gₙₐ = 20
    Eₙₐ = 60
    Vₕₙₐ = -20
    kₙₐ = 15
    gₖ = 10
    Eₖ = -90
    Vₕₖ = -25
    kₖ = 5
    τₙ = 1
    gₘ = 25
    Vₕₘ = -15
    kₘ = 5
    τₕ = 200
    gₙₘ₀ₐ = p[1]
    Eₙₘ₀ₐ = 0
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
           gₙₘ₀ₐ(t) * s∞(V) * (V - Eₙₘ₀ₐ)) / Cₘ + ϵ(process_noise[1])
    du[2] = (n∞(V, Vₕₖ, kₖ) - n) / τₙ + ϵ(process_noise[2])
    du[3] = (h∞(V, Vₕₘ, kₘ) - h) / τₕ + ϵ(process_noise[3])
end

ϵ(process_noise::AbstractFloat) = randn() * process_noise