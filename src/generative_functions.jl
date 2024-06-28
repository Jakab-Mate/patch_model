"""
    create_metabolism(; n_resources::Int64=10, n_levels::Int64=5, energy_yields::String="Uniform_1", seed::Int64=1234)

Generates a universal metabolism

# Optional arguments
- `n_resources::Int64`: Number of possible resources in the system. Default is `10`.
- `n_levels::Int64`: Number of levels of decomposition in the system. Default is `5`.
- `energy_yields::String`: The energy difference between two consecutive levels of decomposition. Default is `1` between all levels. Use "Random" to sample from a uniform distribution between 0 and 2 instead.
- `rng::Int64`: Random number generator seed. Default is `1234`.

# Output
`Stoichiometric matrix, Energy yield matrix`
"""

function create_metabolism(; n_resources::Int64=10, n_levels::Int64=5, energy_yields::String="Uniform_1", seed::Int64=1234)
    rng = MersenneTwister(seed)

    if n_resources < n_levels
        throw(DomainError("Number of resources should be greater than or equal to the number of levels"))
    end

    if n_resources < 2
        throw(DomainError("Number of resources should be greater than 1"))
    end

    levels::Vector{Int64} = vcat(1:n_levels, rand(rng, 2:n_levels-1, n_resources - n_levels))
    levels = sort(levels)
    
    W = zeros(Float64, n_resources)
    for b in n_resources-1:-1:1
        if levels[b] == levels[b+1]
            W[b] = W[b+1]
        else
            if energy_yields == "Uniform_1"
                W[b] = 2 * W[b+1] + 1
            elseif energy_yields == "Random"
                W[b] = 2 * W[b+1] + rand() * 2
            end
        end
    end

    D = zeros(Int64, n_resources, n_resources)
    W_ba = zeros(Float64, n_resources, n_resources)
    for i in 1:n_resources
        for j in 1:n_resources
            if levels[i] > levels[j]
                D[i, j] = 2 ^ (levels[i] - levels[j])
                W_ba[i, j] = W[j] - W[i] * D[i, j]
            end
        end
    end

    return D, W_ba
end

"""
    create_species_pool(D::Matrix; n_families::Int64=5, family_size::Int64=100, dirichlet_hyper::Real=100, between_family_var::Real=0.1, inside_family_var::Real=0.05, h::Real=1, maintenance::Real=0.1, specialist::Real=1, generalist::Real=1, a_dist::Union{Distributions.Sampleable, Nothing}=nothing, k_dist::Union{Distributions.Sampleable, Nothing}=nothing, seed::Int64=1234)

Create a pool of species by sampling reactions from a matrix denoting all possible reactions.

# Mandatory arguments
- `D::Matrix`: Matrix denoting all possible reactions. A square matrix whose size should be equal to the number of possible resources in the system. All reactions will be deemed possible whose values are non-zero.

# Optional arguments
- `n_families::Int64`: Number of families (groups of functionally similar species) in the species pool. Default is `5`.
- `family_size::Int64`: Number of species in each family. Default is `100`.
- `dirichlet_hyper::Real`: Hyperparameter for the Dirichlet distribution that is used for creating the species inside the same family. The higher its value, the more similar they will be. Default is `100`.
- `maintenance::Real`: The expected cost of maintenance accross all species. Default is `0.1`.     
- `between_family_var::Real`: Variance of the normal distribution used to sample the maintenance values between families. Default is `0.1`.
- `inside_family_var::Real`: Variance of the normal distribution used to sample the maintenance values inside families. Default is `0.05`.
- `h::Real`: Controls the allocation of reaction rates inside species. Default is `1`.
- `specialist::Real`: The specialist part of the odds ratio specialists:generalists in the pool. Default is `1`.
- `generalist::Real`: The generalist part of the odds ratio specialists:generalists in the pool. Default is `1`.
- `a_dist::Union{Distributions.Sampleable, Nothing}`: Distribution to sample the strength of host control. Default is `Uniform(0.5, 1.5)`.
- `k_dist::Union{Distributions.Sampleable, Nothing}`: Distribution to sample the critical abundance that triggers host control. Default is `Uniform(99.999, 100.001)`.
- `rng::Int64`: Random number generator seed. Default is `1234`.

# Output
`PoolStruct` with the following fields:
- `pool::Array{Float64, 3}`: The matrices describing the metabolisms of the species inside the species pool.
- `family_ids::Array{Int64}`: The family IDs of species
- `m::Array{Float64}`: The maintenance costs of the species
- `n_reactions::Array{Int64}`: The number of reactions of the species
- `n_splits::Array{Float64}`: Reaction repertoire complexity metric of the species
- `a::Array{Float64}`: The strength of host control on the species
- `k::Array{Float64}`: The critical abundance that triggers host control on the species
"""

function create_species_pool(D::Matrix; n_families::Int64=5, 
    family_size::Int64=100, dirichlet_hyper::Real=100, between_family_var::Real=0.1, inside_family_var::Real=0.05, 
    h::Real=1, maintenance::Real=0.1, specialist::Real=1, generalist::Real=1, 
    a_dist::Union{Distributions.Sampleable, Nothing}=nothing, k_dist::Union{Distributions.Sampleable, Nothing}=nothing, seed::Int64=1234)


    if isnothing(seed)
        rng = MersenneTwister(1234)
    else
        rng = MersenneTwister(seed)
    end

    n_resources::Int64 = length(D[:, 1])
    
    specialist_prob::Int64 = round(Int64, specialist / (generalist + specialist) * 100)
    choices = vcat(fill("Spec", specialist_prob), fill("No", 100-specialist_prob))
    spec_gen = [rand(rng, choices) for _ in 1:n_families]

    pool = zeros(n_resources, n_resources, n_families * family_size)
    family_ids = Array{Int64}(undef, n_families * family_size)
    m = Array{Float64}(undef, n_families * family_size)
    a = Array{Float64}(undef, n_families * family_size)
    k = Array{Float64}(undef, n_families * family_size)
    n_reactions = Array{Int64}(undef, n_families * family_size)
    n_splits = Array{Float64}(undef, n_families * family_size)

    between_family_m_dist = Normal(maintenance, between_family_var)
    if isnothing(a_dist)
        a_dist = Uniform(0.5, 1.5)
    end
    if isnothing(k_dist)
        k_dist = Uniform(99.999, 100.001) # relate this to alpha? (resource input flux)
    end

    for family in 1:n_families
        family_range = (family-1) * family_size + 1 : family * family_size
        family_ids[family_range] .= family
        family_idx_start = (family-1) * family_size 
        family_m_value = rand(rng, between_family_m_dist)
        while family_m_value <= 0
            family_m_value = rand(rng, between_family_m_dist)
        end
        inside_family_m_dist = Normal(family_m_value, inside_family_var)

        family_a_value = rand(rng, a_dist)
        family_k_value = rand(rng, k_dist)
        a[family_range] .= family_a_value
        k[family_range] .= family_k_value

        if spec_gen[family] == "Spec"
            family_number_of_reactions = rand(rng, [4, 5])
        else
            family_number_of_reactions = rand(rng, [2, 3])
        end

        n_reactions[family_range] .= family_number_of_reactions

        prior_values = normalize(rand!(rng, zeros(family_number_of_reactions))) .* dirichlet_hyper
        Dir_dist = Distributions.Dirichlet(prior_values)
        family_reaction_indices = sample_reaction_indices(rng, D, family_number_of_reactions)
        family_n_splits = 0
        for react in family_reaction_indices
            family_n_splits += (react[1] - react[2])
        end
        n_splits[family_range] .= family_n_splits

        for species in 1:family_size
            species_m_value = rand(rng, inside_family_m_dist)
            while species_m_value <= 0
                species_m_value = rand(rng, inside_family_m_dist)
            end
            m[family_idx_start + species] = species_m_value
            species_values = normalize(rand(rng, Dir_dist), h=h)
            for reaction in 1:family_number_of_reactions
                pool[family_reaction_indices[reaction]..., family_size*(family-1)+species] = species_values[reaction]
            end

            #n_splits[family_idx_start + species] = sum((D .- 1) .*  pool[:, :, family_size*(family-1)+species])
        end

    end

    return PoolStruct(pool, family_ids, m, n_reactions, n_splits, a, k)

end