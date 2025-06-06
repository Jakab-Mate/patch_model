# Resource patch model

## A microbial consumer-resource model investigating the microbial communities that form on the surface of different resource patches in the gut

This model builds on our [previous project](https://github.com/Jakab-Mate/MiCroSim.jl), which simulates microbial community dynamics based solely on their metabolic interactions. The model presented here aims to create a more pronounced picture of what happens in the human gut by making a distinction in the underlying chemical reaction system between different type of resource molecules (monomers, simple- and complex resources) using more realistic rules for resource supply (compared to the previous chemostat-like approach), where resources are assumed to appear periodically and in distinct patches due to the feeding of the host. The microbes aggregating on the surface of these patches are modeled as local communities with distinct supply vectors representing the composition of the food particle.

**Please refer to the [project website](https://jakab-mate.github.io/Patches/dev/readme/) for more detailed documentation**

### Installation
To install patch_model directly from the github repository, use:

```julia
using Pkg
Pkg.add(url="https://github.com/Jakab-Mate/patch_model.git")
using Patches
```

### Usage
Detailed instructions for performing simulations can be found in the [User Manual](https://jakab-mate.github.io/MiCroSim.jl/dev/).

### Acknowledgements
While I was the only person working on this code, the model that was implemented mostly consists of the ideas of other researchers. I would like to thank them sincerely for their contributions:

- `Matti Ruuskanen` particularly, for proposing the idea of ephemeral resource patches in the gut.

- `István Scheuring`, `Gergely Boza` and `Leo Lahti` for sharing their valuable knowledge about different factors affecting microbial community structure and metabolism in the gut, along with their ideas about how to represent them in the model.

### Funding

**This project received funding from the European Union’s Horizon 2020 research and innovation programme (under grant agreement 952914; FindingPheno).**

### Contact me

For inquiries and bug reports, contact Jakab Máté: mate.jakab@ecolres.hu


