function sample_pool(p, n_species::Int64=10, n_invaders::Int64=10; seed::Int64=1234)
    """
    Samples species from a species pool

    # Mandatory arguments
    - `p::PoolStruct`: A pool struct containing the pool of species.
    - `n_species::Int64`: Number of species initially present in the habitat. Default is `10`.
    - `n_invaders::Int64`: Number of invading species. Default is `10`.

    # Optional arguments
    - `seed::Int64`: Random number generator seed. Default is `1234`.

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
    rng = MersenneTwister(seed)

    n_resources = size(p.pool, 1)
    n_sampled = n_species + n_invaders

    size_of_pool = size(p.pool, 3)
    if size_of_pool < n_sampled
        throw(DomainError(n_sampled, "\n Number sampled ($n_sampled) should be smaller than the pool ($size_of_pool) \n"))
    end

    if size_of_pool != length(p.family_ids) || size_of_pool != length(p.m) || size_of_pool != length(p.n_reactions) || size_of_pool != length(p.n_splits) || size_of_pool != length(p.a) || size_of_pool != length(p.k)
        throw(DomainError("Size of pool does not match the size of the other pool attributes"))
    end

    species_indices = sample(rng, 1:size_of_pool, n_sampled, replace=false)

    species_C_matrices = p.pool[:, :, species_indices]
    species_family_ids = p.family_ids[species_indices]
    species_m = p.m[species_indices]
    species_n_reactions = p.n_reactions[species_indices]
    species_n_splits = p.n_splits[species_indices]
    species_a = p.a[species_indices] 
    species_k = p.k[species_indices]

    species_initial_abundances = vcat(rand!(rng, zeros(n_species)), zeros(n_invaders))
    resource_initial_abundances = rand!(rng, zeros(n_resources))

    return SampleStruct(n_species, n_invaders, species_C_matrices, species_family_ids, species_m, species_n_reactions, species_n_splits, 
    species_initial_abundances, resource_initial_abundances, species_a, species_k)

end