using Plots
using DynamicalSystems
using DifferentialEquations
using BenchmarkTools
using StatsBase
function loop_lorenz!(du, u, p, t)
    σ, ρ, β = p
    du[1] = σ(t) * (u[2] - u[1])
    du[2] = u[1] * (ρ(t) - u[3]) - u[2]
    du[3] = u[1] * u[2] - β(t) * u[3]
    return nothing

end
transient_T = 50
t = 0:0.01:800
tmax = 800
σ_init = 10.0
β_init = 8 / 3
σ_final = σ_init
β_final = β_init
ρ₀ = 154
ρ₁ = 8
τ = 100

function linear(start::Real, end_::Real, tmax::Real, x::Real)
    return start + (end_ - start) / tmax * x
end

function exponential(start::Real, end_::Real, tmax::Real, x::Real; offset=0.0, τ=0.0)
    if τ == 0
        τ = tmax / log(end_ / start)
    end
    return offset + start * exp(x / τ)
end

params = [t -> linear(σ_init, σ_final, tmax - 600 + transient_T, t - transient_T - 600), t -> exponential(ρ₁, 0, 0, t - transient_T - 600, offset=ρ₀, τ=τ), t -> linear(β_init, β_final, tmax - 600 + transient_T, t - transient_T - 600)]

paper = ContinuousDynamicalSystem(loop_lorenz!, [0.5, 0.5, 0.5], params)
ts = trajectory(paper, 800, Ttr=transient_T)
u = Matrix(ts)
u_std = StatsBase.standardize(ZScoreTransform, u, dims=1)
uS0 = u_std[findall(x -> x <= 600, t), :]
uB0 = u_std[findall(x -> x > 600, t), :]
p = plot3d(uS0[:, 1], uS0[:, 2], uS0[:, 3],
    grid=true,
    lc=:red, label="t\$\\leq\$0", linealpha=1)
plot3d!(p, uB0[:, 1], uB0[:, 2], uB0[:, 3], lc=:black, label="t>0", linealpha=0.8)
display(p)

# ϵ(std::AbstractFloat)::AbstractFloat = randn() * std


# function lorenz_(;u0=[1.0, 0.0, 0.0], ρ=28.0, σ=10.0, β=8 / 3, p=[], t₀=0.0)::ContinuousDynamicalSystem
#     if isempty(p)
#         p = [σ, ρ, β]
#     else
#         p = p
#     end
#     rhs = (du, u, par, t) -> loop_lorenz!(du, u, par, t)
#     return ContinuousDynamicalSystem(rhs, u0, p, t0=t₀)
# end


# condition(u, t, integrator) = true
# affect!(integrator) = set_state!(integrator, get_state(integrator).+ϵ(0.7))
# cb = DiscreteCallback(condition, affect!)


# prob = ODEProblem(loop_lorenz!, [0.5,0.5,0.5],(0,100),[10,28,8/3])
# sol = solve(prob, callback=cb, tstops=[40.0])
# u2 = sol(collect(0:0.01:100))
# p2 = plot3d(u2[1, :], u2[2, :], u2[3, :])
# p2=plot(sol)

# ds = Systems.lorenz()
# diffeq=(alg=Tsit5(), callback=cb)
# u1=[1,0,0]
# u0s=[u1,[0,10,0]]
# pinteg = DynamicalSystems.parallel_integrator(ds,u0s)
# integ=DynamicalSystems.integrator(ds, u1)

    # u2 = trajectory(cde_lorenz, 100.0,diffeq)
# p2 = plot3d(u2[:, 1], u2[:, 2], u2[:, 3])




# function lorenz(du, u, p, t)
#     du[1] = 10.0(u[2] - u[1])
#     du[2] = u[1] * (28.0 - u[3]) - u[2]
#     du[3] = u[1] * u[2] - (8 / 3) * u[3]
# end

# function σ_lorenz(du, u, p, t)
#     du[1] = 3
#     du[2] = 3
#     du[3] = 3
# end

# prob_sde_lorenz = SDEProblem(lorenz, σ_lorenz, [1.0, 0.0, 0.0], (0.0, 10.0))
# sol = solve(prob_sde_lorenz)   
# p1= plot(sol, idxs=(1, 2, 3))

# # u = sol(collect(0:0.01:10))
# p=plot3d(u[1, :], u[2, :], u[3, :],label="sde")
# plot3d!(u2[:, 1], u2[:, 2], u2[:, 3])
# display(p2)