## Reproducibility
The functions in this package (apart from `generic_run()`) are stochastic, and therefore can lead to different results at different times. To ensure reproducibility, the stochastic functions all have a "seed" parameter, which can be used to initialize a specific instance of a random number generator which arrives at the same results every time.

Example:
```julia
using MiCroSim
my_seed = 1234
D, W_ba = create_metabolism(seed=my_seed)
```

## Create a metabolism

```@docs
create_metabolism
```

## Create a species pool

```@docs
create_species_pool
```

## Sample a species pool

```@docs
sample_pool
```

## Simulate dynamics using the sample

```@docs
generic_run
```

## Structs

If you wish to use the structures *PoolStruct* or *SampleStruct*, you will need to specify the module name as well, since these are not exported:
* `MiCroSim.PoolStruct()`
* `MiCroSim.SampleStruct()`