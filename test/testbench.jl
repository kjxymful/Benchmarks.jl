using Benchmarks 
# include("../src/Benchmarks.jl")
# using .Benchmarks

# regimes of bursting neuron
# μ_r = 2.6:0.01:2.7
# μ_s = 9.95:0.001:10.0

# # generate trajectories of bursting neurons with above gₙₘ₀ₐ
# for gₙₘ₀ₐ in μ_r
#     model = bursting_neuron(;gₙₘ₀ₐ)
#     ts =  generate_trajectories(model, 1500, 500; Δt=0.05, process_noise_level=0.0, STD=true, save_name="BN_$(gₙₘ₀ₐ)", PLOT=true, plot_title="BN, gₙₘ₀ₐ=$(gₙₘ₀ₐ)", eval=true, eval_run=0, model_name="BN")
# end

# # generate trajectories of bursting neurons with above gₙₘ₀ₐ in \mu_s
# for gₙₘ₀ₐ in μ_s
#     model = bursting_neuron(;gₙₘ₀ₐ)
#     ts =  generate_trajectories(model, 1500, 500; Δt=0.05, process_noise_level=0.0, STD=true, save_name="BN_$(gₙₘ₀ₐ)", PLOT=true, plot_title="BN, gₙₘ₀ₐ=$(gₙₘ₀ₐ)", eval=true, eval_run=0, model_name="BN")
# end

# # regimes of lorenz
# ρ_r = 23.9:0.01:24

# for ρ in ρ_r
#     model = lorenz(;ρ)
#     ts = generate_trajectories(model, 150, 50; Δt=0.01, process_noise_level=0.0, STD=true, save_name="Lorenz_$(ρ)", PLOT=true, plot_title="Lorenz, ρ=$(ρ)", eval=true, eval_run=0, model_name="Lorenz")
# end

# # regimes of paper lorenz
# ρ_p = 166:0.1:167

# for ρ in ρ_p
#     model = lorenz(;ρ)
#     ts = generate_trajectories(model, 150, 50; Δt=0.01, process_noise_level=0.0, STD=true, save_name="Lorenz_$(ρ)", PLOT=true, plot_title="Lorenz, ρ=$(ρ)", eval=true, eval_run=0, model_name="Lorenz")
# end

t₀ = 0.0
T = 1500.0
transient_T = 500
model_name = "RampUpBN"
model, params = ns_benchmark_systems(model_name, linear, T;t₀, transient_T)
ts = generate_trajectories(model, T, transient_T; Δt=0.05, process_noise_level=0.0, STD=true, save_name="BN_test", PLOT=true, plot_title="BN", eval=true, eval_run=0, model_name=model_name, pl_params =params,TP=true)
