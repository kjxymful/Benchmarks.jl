using DifferentialEquations
using DynamicalSystems
using Benchmarks
using Plots

# include("../src/Benchmarks.jl")
# using .Benchmarks

I=0
Cₘ=6
gₗ=8
Eₗ=-80
gₙₐ=20
Eₙₐ=60
Vₕₙₐ=-20
kₙₐ=15
gₖ=10
Eₖ=-90
Vₕₖ=-25
kₖ=5
τₙ=1
gₘ=25
Vₕₘ=-15
kₘ=5
τₕ=200
gₙₘ₀ₐ=10.1
p = [I, Cₘ, gₗ, Eₗ, gₙₐ, Eₙₐ, Vₕₙₐ, kₙₐ, gₖ, Eₖ, Vₕₖ, kₖ, τₙ, gₘ, Vₕₘ, kₘ, τₕ, gₙₘ₀ₐ]
tspan = (0.0, 10000.0)

prob = ODEProblem(loop_burstn!, [-60.0, 0.0386, 0.0231], tspan, p)
sol = solve(prob)
time = 1000:0.01:10000
u = [sol(t) for t in time]
u = permutedims(hcat(u...))
p1 = plot3d(u[:, 1], u[:, 2], u[:, 3], label="DifferentialEquations")

bn = ContinuousDynamicalSystem(loop_burstn!, [-20,0.5,0.03], p)
u3 = trajectory(bn,10000, Δt=0.01, Ttr=1000)
plot3d!(u3[:,1],u3[:,2],u3[:,3],label="DynamicalSystems")

#display(p1)
# savefig(p1,"test/comparison_integrators.png")

pns = [x->linear(2,4, 1500,x-500)]
ns_loop = ContinuousDynamicalSystem(ns_bursting_neuron!, [-20,0.5,0.03], pns)
u2 = generate_trajectories(ns_loop, 1000,500,PLOT=false)
# u2 = trajectory(ns_loop, 1000, Δt = 0.01, Ttr=500)
p2=plot3d(u2[:,1],u2[:,2],u2[:,3], lc=(:viridis),linez=collect(0:0.01:1000))
# display(p2)

pt = copy(p)
pt[end] = 2
bnt = ContinuousDynamicalSystem(loop_burstn!, [-20,0.5,0.03], pt)
ut = trajectory(bnt, 1500, Δt=0.01,Ttr=5)
ut2 = generate_trajectories(bnt, 1500, 5,PLOT=false)

plot1 = plot3d(ut[:,1],ut[:,2],ut[:,3])
plot2 = plot3d(ut2[:,1],ut2[:,2],ut2[:,3])
plot3 = plot3d!(plot1, ut2[:,1],ut2[:,2],ut2[:,3],label="gen")
end_plot = Plots.plot(plot1, plot2, plot3, layout=3, title=["test" "gen fun" "both"])
# display(end_plot)
@show size(ut2)