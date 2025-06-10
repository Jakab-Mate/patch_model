push!(LOAD_PATH, "../src/")
using Patches, Documenter

ENV["GKSwstype"] = "100"

generated_path = joinpath(@__DIR__, "src")
base_url = "https://github.com/Jakab-Mate/patch_model.jl/blob/main/"
isdir(generated_path) || mkdir(generated_path)

readme_path = joinpath(@__DIR__, "../../main_repo/README.md")

# open(joinpath(generated_path, "readme.md"), "w") do io
#     println(io, """
#     ```@meta
#     EditURL = "$(base_url)README.md"
#     ```
#     """)
#     for line in eachline(readme_path)
#         println(io, line)
#     end
# end

makedocs(
    format = Documenter.HTML(inventory_version = "0.1"),
    authors = "Jakab Máté",
    sitename = "Patches.jl",
    modules = [Patches],
    pages = [
        "Home" => "readme.md",
        "Manual" => "index.md",
        "Tutorials" => Any[
            "Baseline settings" => "example1.md",
            "Order of introductions" => "example2.md"
        ]
    ],
    remotes = Dict(joinpath(pwd(), "main_repo") => Documenter.Remotes.GitHub("Jakab-Mate", "patch_model.jl"))
)

deploydocs(
    repo = "github.com/Jakab-Mate/patch_model.jl",
    branch = "gh-pages",
    devbranch = "main",
    target = "build"
)