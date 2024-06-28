# Add the current project to LOAD_PATH
push!(LOAD_PATH, "@.")

# Load the necessary packages
using MiCroSim, Documenter

# Set up the environment
ENV["GKSwstype"] = "100"

# Generate path and base URL
generated_path = joinpath(@__DIR__, "src")
base_url = "https://github.com/Jakab-Mate/MiCroSim.jl/blob/main/"
isdir(generated_path) || mkdir(generated_path)

# Construct the correct path to README.md
readme_path = joinpath(@__DIR__, "README.md")

# Generate the readme.md file
open(joinpath(generated_path, "readme.md"), "w") do io
    println(io, """
    ```@meta
    EditURL = "$(base_url)README.md"
    ```
    """)
    for line in eachline(readme_path)
        println(io, line)
    end
end

# Generate the documentation
makedocs(
    format = Documenter.HTML(),
    authors = "Jakab Máté",
    sitename = "MiCroSim.jl",
    modules = [MiCroSim],
    pages = [
        "Home" => "readme.md",
        "Manual" => "index.md",
        "Tutorials" => Any[
            "Tutorial 1" => "example1.md",
            "Tutorial 2" => "example2.md"
        ]
    ]
)

# Deploy the documentation
deploydocs(repo = "github.com/Jakab-Mate/MiCroSim.jl", push_preview = true)