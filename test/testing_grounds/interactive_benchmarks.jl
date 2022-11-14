using Benchmarks
using InteractiveDynamics
using DynamicalSystems, GLMakie

# include("src/Benchmarks.jl")
# using .Benchmarks

# sys = bursting_neuron()
# ps = Dict(
#     18 => 2:0.1:12.5
# )
# pnames = Dict(
#     18 => "gₙₘ₀ₐ"
# )
# u1 = [-20,0.1,0.05]
# u2 = [0,0,0]
# u0s = [u1,u2]
# lims = ((-70,0),(0,0.6),(0,0.08))
# interactive_evolution(sys2, u0s; ps, pnames, lims, colors=[:red, :blue], total_span=250)

sys2 = lorenz(ρ=150)
ps = Dict(
    1 => 0:1:30,
    2 => 140:1:180,
    3 => 0:0.5:5
)
pnames = Dict(1 => "σ", 2=>"ρ",3=>"β")
u1=[0.5,0.5,0.5]
u2=[40,10,-10]
u0s=[u1,u2]
interactive_evolution(sys2,u0s;ps,pnames,colors=[:red,:blue],total_span=10)