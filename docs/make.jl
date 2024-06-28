push!(LOAD_PATH, "../src/")
using MiCroSim, Documenter

ENV["GKSwstype"] = "100"

generated_path = joinpath(@__DIR__, "src")
base_url = "https://github.com/Jakab-Mate/MiCroSim.jl/blob/main/"
isdir(generated_path) || mkdir(generated_path)

readme_path = joinpath(@__DIR__, "../../main_repo/README.md")

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

makedocs(
    format = Documenter.HTML(inventory_version = "0.1"),
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
    ],
    remotes = Dict(joinpath(pwd(), "main_repo") => Documenter.Remotes.GitHub("Jakab-Mate", "MiCroSim.jl"))
)

deploydocs(
    repo = "github.com/Jakab-Mate/MiCroSim.jl",
    branch = "gh-pages",
    devbranch = "main",
    target = "build"
)