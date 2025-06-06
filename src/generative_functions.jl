"""
Create a universal metabolism from which species reactions can be sampled

# Arguments
- `n_complex:` The number of different complex (primary) resources
- `n_simple:` The number of different simple resources
- `n_monomer:` The number of different consume_all_monomers
- `limited_pathyways:` Boolean. If true, introduce an additional level of limitation
    (above energetic restraints), where for each resource, some resources can not be produced from it.
    The number of "forbidden" resources is drawn from a uniform distribution
    between 1 and the number of resources that have smaller monomer content than the consumed resource.
- `seed:` Set a specific seed. Default is random.
"""
function create_metabolism(n_complex::Int64,
    n_simple::Int64,
    n_monomer::Int64,
    gapsize::Int64;
    limited_pathways::Bool=false,
    seed::Union{Real, nothing}=nothing)

    if isnothing(seed)
        seed = rand()
    end
    rng = MersenneTwister(seed)
    n_resources = n_complex + n_simple + n_monomer
    monomer_content = zeros(Int64, n_resources)
    for i in 1:n_resources
        if i <= n_complex
            monomer_content[i] = 10 + gapsize
        elseif i <= n_complex + n_simple
            monomer_content[i] = rand(rng, 2:10)
        else
            monomer_content[i] = 1
        end
    end
    sort!(monomer_content, rev=true)

    if limited_pathways
        limiting_matrix = zeros(Int64, n_resources, n_resources)
        for i in 1:n_complex
            allowed_simple_products = sample(rng, n_complex+1:n_complex+n_simple, rand(rng, 1:n_simple), replace=false)
            allowed_monomer_products = sample(rng, n_complex+n_simple+1:n_resources, rand(rng, 1:n_monomer), replace=false)
            all_allowed_products = vcat(allowed_simple_products, allowed_monomer_products)
            for j in all_allowed_products
                limiting_matrix[j, i] = 1
            end
        end

        for i in n_complex+1:n_complex+n_simple
            current_monomer_content = monomer_content[i]
            smaller_monomer_content = []
            for j in n_complex+1:n_complex+n_simple
                if monomer_content[j] < current_monomer_content
                    push!(smaller_monomer_content, j)
                end
            end
            #print("smaller_monomer_content: ", smaller_monomer_content)
            if monomer_content[i] == monomer_content[n_complex + n_simple]
                allowed_simple_products = []
            else
                allowed_simple_products = sample(rng, smaller_monomer_content, rand(rng, 1:length(smaller_monomer_content)), replace=false)
            end

            allowed_monomer_products = sample(rng, n_complex+n_simple+1:n_resources, rand(rng, 1:n_monomer), replace=false)
            all_allowed_products = vcat(allowed_simple_products, allowed_monomer_products)
            for j in all_allowed_products
                limiting_matrix[j, i] = 1
            end
        end
        return monomer_content, limiting_matrix
    else
        return monomer_content, nothing
    end

end


"""
Create a species pool

# Arguments
- `pool_size:` number of species in pool

"""
function create_species_pool(pool_size::Int64, n_complex::Int64, n_simple::Int64, n_monomer::Int64, monomer_content::Array{Int64};
    seed::Union{Real, nothing}=nothing,
    maintenance::Real=0.1, 
    specialist::Real=1, 
    generalist::Real=1, 
    gen_avg_usage::Float64=0.5,
    gen_sd_usage::Float64=0.05,
    spec_avg_usage::Float64=0.5,
    spec_sd_usage::Float64=0.05,
    avg_efficiency::Float64=0.5,
    sd_complex::Float64=0.05,
    sd_simple::Float64=0.1,
    sd_monomer::Float64=0.2,
    gen_n_consumed::Int64=3,
    spec_n_consumed::Int64=2,
    h::Float64=1.0,
    limiting_matrix::Union{Matrix{Int64}, Nothing}=nothing,
    consume_all_monomers::Bool=false)

    # Tradeoff controlling the sum of consumption rates (h)?
    # Linear version of this is encapsulated in the number of reactions tradoff.

    if isnothing(seed)
        seed = rand()
    end
    rng = MersenneTwister(seed)
    n_resources::Int64 = length(monomer_content)
    if isnothing(limiting_matrix)
        limiting_matrix = ones(Int64, n_resources, n_resources)
    end
    resource_sd = vcat(repeat([sd_complex], n_complex), repeat([sd_simple], n_simple), repeat([sd_monomer], n_monomer))
    
    specialist_prob = specialist / (generalist + specialist)
    consumption_rates = []
    energy = []
    production = []

    types = Array{String}(undef, pool_size)
    m = Array{Float64}(undef, pool_size)
    n_reactions = Array{Int64}(undef, pool_size)
    n_splits = Array{Float64}(undef, pool_size)
    
    println("generating species")
    for species in 1:pool_size
        decider = rand(rng)
        if consume_all_monomers
            if decider < specialist_prob
                types[species] = "spec"
                n_reactions[species] = spec_n_consumed
                consumed_indices = sample(rng, n_complex+1:n_complex+n_simple, spec_n_consumed, replace=false)
                species_consumption_rates, species_energy, species_production, species_n_splits = partition_resources(
                    n_resources, consumed_indices, 1/spec_n_consumed, monomer_content, limiting_matrix, spec_avg_usage, spec_sd_usage, avg_efficiency, resource_sd, rng)
                species_consumption_rates[n_complex+n_simple+1:n_resources] .= 1.0
            else
                types[species] = "gen"
                n_reactions[species] = gen_n_consumed
                consumed_primary = sample(rng, 1:n_complex, 1)
                primary_and_simple_range = [x for x in 1:n_complex+n_simple if x != consumed_primary[1]]
                consumed_simple = sample(rng, primary_and_simple_range, gen_n_consumed-1, replace=false)
                consumed_indices = vcat(consumed_primary, consumed_simple)
                species_consumption_rates, species_energy, species_production, species_n_splits = partition_resources(
                    n_resources, consumed_indices, 1/gen_n_consumed, monomer_content, limiting_matrix, gen_avg_usage, gen_sd_usage, avg_efficiency, resource_sd, rng)
                species_consumption_rates[n_complex+n_simple+1:n_resources] .= 1.0
            end
        else
            if decider < specialist_prob
                types[species] = "spec"
                n_reactions[species] = spec_n_consumed
                consumed_indices = sample(rng, n_complex+1:n_resources, spec_n_consumed, replace=false)
                species_consumption_rates, species_energy, species_production, species_n_splits = partition_resources(
                    n_resources, consumed_indices, 1/spec_n_consumed, monomer_content, limiting_matrix, spec_avg_usage, spec_sd_usage, avg_efficiency, resource_sd, rng)
            else
                types[species] = "gen"
                n_reactions[species] = gen_n_consumed
                consumed_primary = sample(rng, 1:n_complex, 1)
                consumed_simple = sample(rng, n_complex+1:n_complex+n_simple+n_monomer, gen_n_consumed-1, replace=false)
                consumed_indices = vcat(consumed_primary, consumed_simple)
                species_consumption_rates, species_energy, species_production, species_n_splits = partition_resources(
                    n_resources, consumed_indices, 1/gen_n_consumed, monomer_content, limiting_matrix, gen_avg_usage, gen_sd_usage, avg_efficiency, resource_sd, rng)
            end
        end

        if consume_all_monomers
            #println("species_consumption_rates before normalization: ", species_consumption_rates)
            species_consumption_rates[1:n_complex+n_simple] = normalize(species_consumption_rates[1:n_complex+n_simple], h=h)
            #println("species_consumption_rates after normalization: ", species_consumption_rates)
        else   
            #println("species_consumption_rates before normalization: ", species_consumption_rates)
            species_consumption_rates = normalize(species_consumption_rates, h=h)
            #println("species_consumption_rates after normalization: ", species_consumption_rates)
        end

        species_m_value = bounded_rand(rng, Normal(maintenance, 0.01), (0, 1))
        m[species] = species_m_value

        push!(consumption_rates, species_consumption_rates)
        push!(energy, species_energy)
        push!(production, species_production)
        n_splits[species] = species_n_splits
    end

    println("species generated")

    return PoolStruct(consumption_rates, energy, production, m, n_reactions, n_splits, types)
end