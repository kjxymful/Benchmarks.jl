function interactive_benchmarks(args::Dict{String,Any})
    sys_name = args["name"]
    dynamical_noise = args["interactive_noise"]
    if occursin("Lorenz", sys_name) || occursin("lorenz",sys_name)
        u1=[0.5,0.5,0.5]
        u2=[40,10,-10]
        u0s=[u1,u2]
        lims = ((-20,40),(-25,25),(-10,60))
        total_span=10f0
        if occursin("Paper", sys_name)
            lims=((-30,45),(-120,80),(0,250))
        end
        sys = lorenz_pn()
    elseif occursin("BN",sys_name) || occursin("Bursting", sys_name) || occursin("bursting",sys_name)
        sys = bursting_neuron_pn()
        u1 = [-20,0.1,0.05]
        u2 = [-21,0.1,0.05]
        u0s = [u1,u2]
        lims = ((-70,0),(0,0.6),(0,0.08))
        total_span=250f0
    else
        throw("Not Implemented $sys_name")
    end
    std_ = dynamical_noise ? std_model(sys,total_span,Int(total_span÷10)) : [0,0,0] 
    ps, pnames = choose_p_dict(sys_name, dynamical_noise)
    interactive_evolution(sys, u0s; ps, pnames, lims, colors=[:red, :blue], total_span,diffeq=(alg=Tsit5(), callback=dynamical_noise_parallel_callback(dynamical_noise,std_)))
end

function choose_p_dict(sys_name::String,dynamical_noise::Bool)
    if occursin("Lorenz", sys_name) || occursin("lorenz",sys_name)
        ps = Dict(
            1 => 0:1:30,
            2 => 0:1:40,
            3 => 0:0.1:5
        )
        pnames = Dict(1 => "σ", 2=>"ρ",3=>"β")
        if occursin("Paper", sys_name)
            ps = Dict(
                1 => 0:1:30,
                2 => 140:1:180,
                3 => 0:0.1:5
            )
        end  
        if dynamical_noise
            push!(ps, (4=>0:0.05:0.5))
            push!(pnames,(4=>"DNLvl"))
        end
    elseif occursin("BN",sys_name) || occursin("Bursting", sys_name) || occursin("bursting",sys_name)
        ps = Dict(
            18 => 2:0.1:12.5
        )
        pnames = Dict(
            18 => "gₙₘ₀ₐ"
        )
        if dynamical_noise
            push!(ps, (19=>0:0.05:0.5))
            push!(pnames,(19=>"DNLvl"))
        end
    else
        throw("$System not Implemented")
    end
    return ps, pnames
end