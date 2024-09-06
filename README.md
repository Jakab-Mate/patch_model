# Microbial Cross-feeding Community Simulator

[![Build Status](https://github.com/Jakab-Mate/MiCroSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Jakab-Mate/MiCroSim.jl/actions/workflows/CI.yml?query=branch%3Amain)

## A general-purpose microbial consumer-resource model that outputs simulated data in the SummarizedExperiment (SE) format.

This project offers a quick and easy method for simulating microbial community dynamics based solely on their metabolic interactions. The model uses matrix representations for species, which encode the metabolites a species can consume, as well as those it can produce. Community dynamics are simulated using a set of Ordinary Differential Equations (ODEs), and the produced time-series data is stored in a SummarizedExperiment data container.

**Please refer to the [project website](https://jakab-mate.github.io/MiCroSim.jl/dev/readme/) for more detailed documentation**

### Applicability

* **Simulate the growth of an initially populated microbial community and find the equilibrium abundances**

* **Simulate the succession of an initially empty microbial habitat where colonizers appear one by one**

* **Simulate the resilience of a microbial community against invaders appearing one by one**

### Installation
To install MiCroSim.jl directly from the github repository, use:

```julia
using Pkg
Pkg.add(url="https://github.com/Jakab-Mate/MiCroSim.jl.git")
```

### Workflow
The functions in this package rely on each other's outputs, so generally, you will want to use them in the following order:
1. **create_metabolism(...)**<span style="display:inline-block; width: 56px;"></span>Generates the set of possible reactions (net conversions)

2. **create_species_pool(...)**<span style="display:inline-block; width: 48px;"></span>Generates the pool of possible species by sampling from the possible reactions

3. **sample_pool(...)**<span style="display:inline-block; width: 100px;"></span>Chooses the species that will appear in the simulation by sampling from the species pool

4. **generic_run(...)**<span style="display:inline-block; width: 107px;"></span>Takes the sampled species and simulates their dynamics based on a set of ODEs

4. **spatial_run(...)**<span style="display:inline-block; width: 107px;"></span>Spatial alternative of `generic_run`, models multiple interconnected local communities

For more detailed explanations about what each function does, please refer to the [workflow section of the home page](https://jakab-mate.github.io/MiCroSim.jl/dev/readme/#Workflow)

Detailed instructions for using each function can be found in the [User Manual](https://jakab-mate.github.io/MiCroSim.jl/dev/).

### Design your own metabolism

The best way to contribute to this project is by curating universal metabolisms in the form of stoichiometric and energy yield matrices. Admittedly, the reaction systems that may arise from `create_metabolism()` are limited, but more complex metabolic networks can also be implemented, for example modeling synthetic processes by setting energy yields negative (that is, a species invests into producing a metabolite). Furthermore, pathway databases such as KEGG coupled with microbial whole genome data open the possibility for deriving net conversions from real-world experiments.

### Acknowledgements
This project has benefited from contributions and insights of the following individuals and groups:

* **István Scheuring** and **Gergely Boza** from the Centre for Ecological Research, Budapest, provided essential theoretical guidance for the model's development.

* **Leo Lahti**, **Giulio Benedetti**, **Moien Khalighi** and the rest of the [**TurkuDataScience**](https://datascience.utu.fi/) team form the University of Turku, Turku, were instrumental in setting up and optimizing the Julia package. Perhaps more importantly, I would like to thank them for their hospitality throughout my secondment in Turku.

* The model was inspired by the work of [Goldford et al. (2018)](https://doi.org/10.1126/science.aat1168).

### Funding

**This project received funding from the European Union’s Horizon 2020 research and innovation programme (under grant agreement 952914; FindingPheno).**

### Contact me

For inquiries and bug reports, contact Jakab Máté: mate.jakab@ecolres.hu


