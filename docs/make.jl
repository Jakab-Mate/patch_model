push!(LOAD_PATH, "../src/")
using Pkg

Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Patches
using Documenter

ENV["GKSwstype"] = "100"

makedocs(
    format = Documenter.HTML(inventory_version = "0.1"),
    authors = "Jakab Máté",
    sitename = "patch_model",
    pages = [
        "Home" => "readme.md",
        "Manual" => "index.md",
        "Tutorials" => Any[
            "Default parameters" => "example1.md",
            "Order of introductions" => "example2.md"
        ]
    ],
    repo = Documenter.Remotes.GitHub("Jakab-Mate", "patch_model")
)

deploydocs(
    repo = "github.com/Jakab-Mate/patch_model.git",
    branch = "gh-pages",
    devbranch = "main",
    target = "build",
    deploy_config = Documenter.GitHubActions()
)