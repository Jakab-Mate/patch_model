module Patches

using Random: rand!, MersenneTwister, shuffle
using Dates: now
using Plots #: plot, plot!, savefig
theme(:wong)
using StatsBase: sample
using MicrobiomeAnalysis
using Distributions: Normal, Uniform, Dirichlet
using DifferentialEquations: ODEProblem, solve, KenCarp4, CallbackSet, PeriodicCallback
using Graphs: SimpleDiGraph, add_edge!, add_vertex!, vertices, edges
using GraphPlot: gplot, circular_layout
using Compose: compose, draw, rectangle, fill, cm
using Base.Threads
using JLD2: @save
import Cairo, Fontconfig

include("structs.jl")
include("generative_functions.jl")
include("sample_pool.jl")
include("equations.jl")
include("callbacks.jl")
include("run_functions.jl")
include("plot_functions.jl")
include("helper_functions.jl")
include("spatial_equations.jl")
include("essential_R_equations.jl")
include("analyze_ts.jl")
include("repeat_params.jl")
include("grab_ends.jl")

export(generic_run)
export(spatial_run)
export(sample_pool)
export(create_metabolism)
export(create_species_pool)
export(repeat_params)
export(grab_ends)

end
