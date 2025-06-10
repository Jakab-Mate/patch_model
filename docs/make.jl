push!(LOAD_PATH, "../src/")
using Patches, Documenter

ENV["GKSwstype"] = "100"

generated_path = joinpath(@__DIR__, "src")
base_url = "https://github.com/Jakab-Mate/patch_model/blob/main/"
isdir(generated_path) || mkdir(generated_path)

readme_path = joinpath(@__DIR__, "../../main_repo/README.md")

makedocs(
    format = Documenter.HTML(inventory_version = "0.1"),
    authors = "Jakab Máté",
    sitename = "patch_model",
    modules = [Patches],
    pages = [
        "Home" => "readme.md",
        "Manual" => "index.md",
        "Tutorials" => Any[
            "Baseline settings" => "example1.md",
            "Order of introductions" => "example2.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/Jakab-Mate/patch_model",
    branch = "gh-pages",
    devbranch = "main",
    target = "build"
)