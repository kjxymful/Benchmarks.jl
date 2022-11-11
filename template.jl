using Benchmarks 
using DynamicalSystems
using Plots


T = 150
transient_T = 10
u0 = [0.5,0.5,0.5] 
# Δt is another kwarg

BN = bursting_neuron(gₙₘ₀ₐ=9, process_noise=zeros(3), u0=u0) # same for lorenz
ts = generate_trajectories(BN, T, transient_T,PLOT=false)

# non-stationary 
jacobian = ns_lorenz!

σ(t) = linear(10,15,T,t-transient_T)
ρ(t) = linear(20,20,T,t-transient_T)
β(t) = linear(8/3,8/3,T,t-transient_T)
parameters = [σ,ρ,β]

ns_Lorenz = ns_systems(jacobian,parameters, u0, "temlate")
ts = generate_trajectories(ns_Lorenz, T, transient_T,PLOT=false)
# or 
ns_CLorenz = ContinuousDynamicalSystem(jacobian, u0,parameters,t0=0.0)
ts = generate_trajectories(ns_CLorenz, T, transient_T,PLOT=false)
# or
ts = trajectory(ns_CLorenz, T, Ttr=transient_T)
