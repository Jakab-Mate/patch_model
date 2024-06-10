module TestPkg

using Random
using Dates
using Plots 
using Distributions
using DifferentialEquations 
using LinearAlgebra 
using MicrobiomeAnalysis

include("structs.jl")
include("generative_functions.jl")
include("sample_pool.jl")
include("equations.jl")
include("callbacks.jl")
include("run_functions.jl")
include("plot_functions.jl")
include("helper_functions.jl")

export(generic_run)
export(sample_pool)
export(create_metabolism)
export(create_species_pool)

end
