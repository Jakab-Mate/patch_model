### Simulate spatial dynamics

```julia
using MiCroSim

# Create some metabolism
D, W_ba = create_metabolism(n_resources=n_resources)

# Create a species pool
pool = create_species_pool(D, n_families=5, family_size=10)

# Draw a sample of 10 species and 10 invaders
sample = sample_pool(pool, 10,  10)

# Simulate dynamics of a network that is partitioned into 3 local communities
out = spatial_run(3, sample, D=D, W_ba=W_ba, phi=0.1, eta=0.1, t_span=(0, 400), host_regulation=false)
```