using Benchmarks

# see doc for available options

function main(args::Dict{String,Any})
    if args["interactive"]
        interactive_benchmarks(args)
    else
        generate_benchmarks(args)
    end
end

main(parse_commandline())
