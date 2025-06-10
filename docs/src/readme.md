# Resource patch model

This model builds on our [previous project](https://github.com/Jakab-Mate/MiCroSim), which simulates microbial community dynamics based solely on their metabolic interactions. The model presented here aims to create a more pronounced picture of what happens in the human gut by making a distinction in the underlying chemical reaction system between different type of resource molecules (monomers, simple- and complex resources) using more realistic rules for resource supply (compared to the previous chemostat-like approach), where resources are assumed to appear periodically and in distinct patches due to the feeding of the host. The microbes aggregating on the surface of these patches are modeled as local communities with distinct supply vectors representing the composition of the food particle.

Moreover, we distinguish a mucus community, motivated by the fact that the mucus operates very differently from the ephemeral resource patches of the lumen: The mucus resource is less energy-rich than other complex resources, but it provides a steady supply. Therefore, the mucus can serve as a refuge in times of starvation. We think it is more realistic, to have the mucus resource be made out of different building blocks than the rest of the complex resources, but a switch parameter of the model (`mucus_disjunct`) allows for either option.

Local communities, including the mucus, follow the same dynamics: species can consume resources, and for each resource consumed, some amount of its monomer content contributes to population growth, while the rest contributes to byproduct formation. Byproducts are released into a shared environment, where they are available for consumption by other species. Local communities are connected by diffusion, and (optionally) advection.

This model is one of succession: we periodically introduce invaders a pre-defined number of times, which can appear in either local community. When invasions are fast, the system is in a constant transient state, therefore we need some additional time to see the equilibrium state community after invasion events stop. Several measures and plots are calculated form the generated time-series data, see the Manual for details.

## Acknowledgements

While I was the only person working on this code, the model that was implemented mostly consists of the ideas of other researchers. I would like to thank them sincerely for their contributions:

* Matti Ruuskanen particularly, for proposing the idea of ephemeral resource patches in the gut.

* István Scheuring, Gergely Boza and Leo Lahti for sharing their valuable knowledge about different factors affecting microbial community structure and metabolism in the gut, along with their ideas about how to represent them in the model.

## Funding

This project received funding from the European Union’s Horizon 2020 research and innovation programme (under grant agreement 952914; FindingPheno).

## Contact me

For inquiries and bug reports, contact Jakab Máté: mate.jakab@ecolres.hu


