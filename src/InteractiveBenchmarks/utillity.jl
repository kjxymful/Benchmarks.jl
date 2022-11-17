using Statistics: std
function std_model(ds::ContinuousDynamicalSystem,tmax, transient_T::Int; Δt=0.01)
    tseries = Matrix(trajectory(ds, tmax, Ttr=transient_T,Δt=Δt))
    std_ = [std(tseries[:, i]) for i in axes(tseries, 2)]
    return std_
end

function dynamical_noise_parallel_callback(dynamical_noise::Bool, std::AbstractVector)::DiscreteCallback
    condition(u, t, integrator) = dynamical_noise
    function affect!(integrator)
        noise_level = integrator.p[end].*std
        for i in 1:2
            set_state!(integrator, get_state(integrator, i) .+ ϵ.(noise_level), i)
        end
    end
    cb = DiscreteCallback(condition, affect!)
    return cb
end

ϵ(dynamical_noise_level::AbstractFloat)::AbstractFloat = randn() * dynamical_noise_level
