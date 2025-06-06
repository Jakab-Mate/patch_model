# Resource patch model

[![Build Status](https://github.com/Jakab-Mate/MiCroSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Jakab-Mate/MiCroSim.jl/actions/workflows/CI.yml?query=branch%3Amain)

## A microbial consumer-resource model investigating the microbial communities that form on the surface of different resource patches in the gut

This project offers a quick and easy method for simulating microbial community dynamics based solely on their metabolic interactions. The model uses matrix representations for species, which encode the metabolites a species can consume, as well as those it can produce. Community dynamics are simulated using a set of Ordinary Differential Equations (ODEs), and the produced time-series data is stored in a SummarizedExperiment data container.

**Please refer to the [project website](https://jakab-mate.github.io/MiCroSim.jl/dev/readme/) for more detailed documentation**

### Applicability

* **Simulate the growth of an initially populated microbial community and find the equilibrium abundances**

### Installation
To install patch_model directly from the github repository, use:

```julia
using Pkg
Pkg.add(url="https://github.com/Jakab-Mate/patch_model.git")
using Patches
```

### Workflow

For more detailed explanations about what each function does, please refer to the [workflow section of the home page](https://jakab-mate.github.io/MiCroSim.jl/dev/readme/#Workflow)

Detailed instructions for using each function can be found in the [User Manual](https://jakab-mate.github.io/MiCroSim.jl/dev/).

### Acknowledgements

### Funding

**This project received funding from the European Union’s Horizon 2020 research and innovation programme (under grant agreement 952914; FindingPheno).**

### Contact me

For inquiries and bug reports, contact Jakab Máté: mate.jakab@ecolres.hu


