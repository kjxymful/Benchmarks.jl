function gen_series(ds::GeneralizedDynamicalSystem, num_T::Int, ΔT::AbstractFloat, transient_T::Int)
    T_end = num_T*ΔT
    ts = trajectory(ds, T_end; Δt=ΔT, Ttr=transient_T)
    return Matrix(ts)
end
