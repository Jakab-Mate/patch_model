module Patches

using Random: rand!, MersenneTwister, shuffle
using Dates: now
using Plots #: plot, plot!, savefig
theme(:wong)
using LinearAlgebra: norm
using StatsBase: sample
using MicrobiomeAnalysis
using Distributions: Normal, Exponential, Uniform, Beta, Gamma, Dirichlet
using DifferentialEquations: ODEProblem, solve, Tsit5, Rodas4
using Graphs: Graph, SimpleGraph, DiGraph, SimpleDiGraph, add_edge!, add_vertex!, vertices, edges
using GraphPlot: graphplot, plot, plot!, savefig
using Compose: Compose, @draw, @cairo, @cairo_dpi, @cairo_dpi_size
using CSV: CSV, read, write
using Base.Threads
using JLD2: jldopen, JLD2, save, load, @load, @save
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
