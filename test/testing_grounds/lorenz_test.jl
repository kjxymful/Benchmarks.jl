using DynamicalSystems
using Benchmarks
using Plots

σ_final = 5.0
ρ_final = 23.0
β_final = 0.5

# include("../src/Benchmarks.jl")
# using .Benchmarks

σ = x->linear(10f0, σ_final, 150f0,x-50)
ρ = x->linear(28f0, ρ_final, 150f0, x-50)
β = x->linear(8/3, β_final, 150f0, x-50)

p = [σ,ρ,β]
sl = ContinuousDynamicalSystem(ns_lorenz!, [0.5 for i=1:3], p, t0=-100.0)
ts = trajectory(sl, 150, Δt = 0.01, Ttr=50)
@show size(ts)
time = collect(0:0.01:150)
p_st = plot3d(ts[:,1],ts[:,2],ts[:,3], lc=(:viridis), linez=collect(0:0.01:150))
display(p_st)
# p_par = plot(time, σ.(time), label="σ")
# plot!(time,ρ.(time), label="ρ")
# plot!(time, β.(time), label="β")
# l = grid(2,1, heights=[0.8,0.2])
# pl = plot(p_st, p_par, layout=l)
# display(pl)
