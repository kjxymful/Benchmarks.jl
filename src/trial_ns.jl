
function trial_ns_lorenz(num_trials::Int, trial_length::Int, Δt::AbstractFloat, μ_start::Float32, μ_end::Float32;process_noise_level=0.0,transient_T=0)
    ρ_diff = (μ_end - μ_start)/num_trials
    ρs = μ_start:ρ_diff:μ_end

    tmax = trial_length*Δt
    u0 = [-10.0,20.0,-15.0]
    trial_lorenz = Vector{AbstractMatrix}()
    for ρ in ρs
        ds = lorenz(u0=u0, ρ=ρ)
        tseries = generate_trajectories(ds, tmax, transient_T; Δt,process_noise_level, PLOT=false)
        u0 = tseries[end,:]
        transient_T = 0
        tmax = trial_length*Δt
        push!(trial_lorenz, tseries)
    end
    return trial_lorenz
end

function split_lorenz_trials(System::String, num_trials::Int, seq_length::Int, transient_T::Int, Δt::AbstractFloat; process_noise_level=0.0)
    num_T = num_trials*seq_length
    tmax = num_T * Δt

    ns_model,params = ns_benchmark_systems(System, linear, tmax)

    tseries = generate_trajectories(ns_model, tmax, transient_T; Δt, process_noise_level, model_name=System, PLOT=false)

    trial_lorenz = Vector{AbstractMatrix}()

    for i in 0:(num_trials-1)
        push!(trial_lorenz, tseries[i*seq_length+1:i*seq_length+seq_length+1,:])
    end
    return trial_lorenz
end
    