using Random, Dates, Plots, DifferentialEquations, LSODA, LinearAlgebra, MicrobiomeAnalysis

include("structs.jl")
include("generative_functions.jl")
include("sample_pool.jl")
include("equations.jl")
include("callbacks.jl")
include("run_functions.jl")
include("plot_functions.jl")

# seed = now().instant.periods.value
# rng = MersenneTwister(seed)

# n_resources = 10
# n_species = 5
# n_invaders = 100

# D, W_ba = generative_functions.create_metabolism(rng, n_resources=n_resources, energy_yields="Random")
# pool = generative_functions.create_species_pool(rng, D, n_families=20)

# s = sample_pool(rng, pool, n_resources, n_species, n_invaders)

# out = generic_run("/home/jakab/Julia/CR_results", D, W_ba, s,
#     n_resources=n_resources,
#     n_species=n_species,
#     n_invaders=n_invaders,
#     t_span=(0, 3000))

using PkgTemplates
t = Template(user="Jakab-Mate")
t("TestPkg")
