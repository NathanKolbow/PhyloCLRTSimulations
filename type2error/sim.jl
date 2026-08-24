length(ARGS) > 0 || error("Must provide 1 command line argument for the value of ϵ (0, 0.0001, 0.001, 0.01, 0.1, 0.5, or 1.0)")
ϵ = parse(Float64, ARGS[1])
ϵ in [0, 0.0001, 0.001, 0.01, 0.1, 0.5, 1.0] || error("ϵ=$(ϵ) not allowed.")

include("../includes.jl")
