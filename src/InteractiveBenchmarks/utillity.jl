
function std_model(ds::ContinuousDynamicalSystem,tmax, transient_T::Int; Δt=0.01, t₀=0.0)
    tseries = Matrix(trajectory(ds, tmax, Ttr=transient_T,Δt=Δt, t0=t₀))
    std_ = [std(tseries[:, i]) for i in axes(tseries, 2)]
    return std_
end

function dynamical_noise_parallel_callback(dynamical_noise_level::AbstractFloat, std::AbstractVector)::DiscreteCallback
    noise_level = dynamical_noise_level .* std
    dn = dynamical_noise_level == 0 ? false : true
    condition(u, t, integrator) = dn
    function affect_pint!(integrator)
        for i in 1:2
            set_state!(integrator, get_state(integrator, i) .+ ϵ.(noise_level), i)
        end
    end
    cb = DiscreteCallback(condition, affect_pint!)
    return cb
end