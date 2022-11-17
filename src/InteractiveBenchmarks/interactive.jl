function interactive_benchmarks(args::Dict{String,Any})
    sys_name = args["name"]
    pn = args["interactive_noise"]
    if occursin("Lorenz", sys_name) || occursin("lorenz",sys_name)
        ps = Dict(
            1 => 0:1:30,
            2 => 0:1:40,
            3 => 0:0.1:5
        )
        pnames = Dict(1 => "σ", 2=>"ρ",3=>"β")
        u1=[0.5,0.5,0.5]
        u2=[40,10,-10]
        u0s=[u1,u2]
        lims = ((-20,40),(-25,25),(-10,60))
        total_span=10f0
        if occursin("Paper", sys_name)
            std = [20,32,30]
            ps = Dict(
                1 => 0:1:30,
                2 => 140:1:180,
                3 => 0:0.1:5
            )
            lims=((-30,45),(-120,80),(0,250))
        end
        sys = lorenz()
    elseif occursin("BN",sys_name) || occursin("Bursting", sys_name) || occursin("bursting",sys_name)
        sys = bursting_neuron()
        ps = Dict(
            18 => 2:0.1:12.5
        )
        pnames = Dict(
            18 => "gₙₘ₀ₐ"
        )
        u1 = [-20,0.1,0.05]
        u2 = [0,0,0]
        u0s = [u1,u2]
        lims = ((-70,0),(0,0.6),(0,0.08))
        total_span=250f0
    else
        throw("Not Implemented $sys_name")
    end
    std_ = pn==0 ? [0,0,0] : std_model(sys,total_span,Int(total_span÷10))
    interactive_evolution(sys, u0s; ps, pnames, lims, colors=[:red, :blue], total_span,diffeq=(alg=Tsit5(), callback=dynamical_noise_parallel_callback(pn,std_)))
end

