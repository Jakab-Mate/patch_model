"""
    sample_pool(p::PoolStruct, n_species::Int64, n_invaders::Int64; seed::Int64=1234)

Samples species from a species pool

# Mandatory arguments
- `p::PoolStruct`: A pool struct containing the pool of species.
- `n_species::Int64`: Number of species initially present in the habitat.
- `n_invaders::Int64`: Number of invading species.

# Optional arguments
- `seed::Int64`: Random number generator seed. Default is `1234`.
- `n_comms::Int64`: Number of communities. Default is `1`.
- `ph::Bool`: Whether to sample pH values. If true, optimal pH value is sampled from the (0, 14) range for each species. 
Default is `false`, where optimal pH is set to 7.0 for all species.

# Output
`SampleStruct` with the following fields:
- `n_species::Int64`: Number of species initially present in the habitat.
- `n_invaders::Int64`: Number of invading species.
- `C::Array{Float64, 3}`: The matrices describing the metabolism of the sampled species.
- `family_ids::Array{Int64}`: The family IDs of the sampled species.
- `m::Array{Float64}`: The maintenance costs of the sampled species.
- `n_reactions::Array{Int64}`: The number of reactions of the sampled species.
- `n_splits::Array{Float64}`: Reaction repertoire complexity metric of the sampled species.
- `species_abundance::Array{Float64}`: The initial abundances of the sampled species.
- `resource_abundance::Array{Float64}`: The initial abundances of resources.
- `a::Array{Float64}`: The strength of host control on the sampled species.
- `k::Array{Float64}`: The critical abundance that triggers host control on the sampled species.
"""
function sample_pool(p::PoolStruct, n_species::Int64, n_invaders::Int64; n_comms::Int64=1, seed::Int64=1234, ph::Bool=false, fix_order::Bool=false)
    rng = MersenneTwister(seed)
    n_resources = length(p.consumption_rates[1])
    n_sampled = n_species + n_invaders
    size_of_pool = length(p.consumption_rates)
    if size_of_pool < n_sampled
        throw(DomainError(n_sampled, "\n Number sampled ($n_sampled) should be smaller than the pool ($size_of_pool) \n"))
    end

    species_indices = sample(rng, 1:size_of_pool, n_sampled, replace=false)
    if !fix_order
        species_indices = shuffle(species_indices)
    end
    species_consumption_rates = p.consumption_rates[species_indices]
    species_energy = p.energy[species_indices]
    species_production_matrices = p.production_matrices[species_indices]
    species_n_reactions = p.n_reactions[species_indices]
    species_n_splits = p.n_splits[species_indices]
    species_m = p.m[species_indices]
    species_types = p.types[species_indices]

    species_initial_abundances = vcat(rand!(rng, zeros(n_species)), zeros(n_invaders))
    resource_initial_abundances = rand!(rng, zeros(n_resources))

    if ph
        species_ph_opts = rand!(rng, zeros(n_sampled)) .* 14.0
    else
        species_ph_opts = fill(7.0, n_sampled)
    end

    # All communities have random initial abundances (?)
    if n_comms > 1
        for i in 2:n_comms
            species_initial_abundances = vcat(species_initial_abundances, zeros(n_species+n_invaders))
            resource_initial_abundances = vcat(resource_initial_abundances, zeros(n_resources))
        end
    end

    println("Length of resource states: ", length(resource_initial_abundances))

    return SampleStruct(n_species, n_invaders, species_consumption_rates, species_energy, species_production_matrices, species_m, species_n_reactions,
        species_n_splits, species_initial_abundances, resource_initial_abundances, species_ph_opts, species_types)

end