
To run simulations, we resort to a single function with a wide variety of parameters. 

```@docs
Patches.repeat_params
```

Each simulation will have its own subfolder in the specified output folder, where the following things are saved:

* `params`: Parameters of the simulation
* `beta_diversity plot`: A measure of how different the composition of local communities are
* `shannon_diversity plot`: A measure of the diversity of all populations of the simulation pooled together
* `MCI plot`: Metabolic Capacity Index of the of the system. Tracks the total number of different reactions by taking all of the species stoichiometric matrices and adding 1 for each position where at least 1 species has a non-zero entry
* `consumer_types plot`: Tracks the total number of primary consumer (generalist) and secondary consumer (specialist) species across the system
* `species_counts plot`: Tracks the total number of different species across the system
* `comm_totals plot`: Tracks the total abundance of members for each local community
* `community graph plots`: A separate graph for each local community, depicting the active pathways. Species nodes are blue, resource nodes are green. Red arrows indicate consumption, green arrows indicate production. Arrow widths are weighted according to the pathway fluxes.

Outside of the individual run folders, we will find the following outputs:

* `Aggregated versions` of the Beta diversity, Shannon diversity, Species counts, MCI, Consumer types and Community totals figures. Standard deviation is shown in all but the last two figures of the above list. 

* `Tr_rel_something plots` correspond to the same measures, but they take an average of these measures around the end of the transient (tr) and at the end of the relaxation (rel) period for each simulation. The mean and standard deviation are shown.

* `JLD2 format files` are used to store the data needed for generating all of the above plots.