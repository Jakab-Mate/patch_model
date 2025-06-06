
"""
Return the indices of non-zero elements in a matrix as a tuple of tuples.
"""
function non_zero_pos(A)
    indices = findall(x -> x != 0, A)
    tuple_indices = Tuple.(indices)
    return tuple_indices
end

"""
Normalize an array `A` by dividing each element by the sum of the elements.
If `h` is provided, normalize by the sum of the elements raised to the power of `h`.
"""
function normalize(A; h::Real=0)
    if isempty(A)
        throw(DomainError("Cannot normalize an empty array."))
    else
        if sum(A) == 0
            throw(DomainError("Cannot normalize an array that sums to zero."))
        end
    end

    if h == 0
        return A ./ sum(A)
    else
        return A ./ sum(A .^ h)
    end     
end

"""
Sample `number_of_reactions` amount of positions from a `D` metabolism matrix.
The function will sample from the non-zero positions of `D` and return the indices of the sampled reactions.
"""
function sample_reaction_indices(rng, D, number_of_reactions)
    # TODO QUESTION: Allow simultaneous consumtion and production of the same resource?
    total_different_reactions = length(non_zero_pos(D))
    if total_different_reactions < number_of_reactions
        throw(DomainError("Number of sampled reactions ($number_of_reactions) should be smaller than the total number of reactions ($total_different_reactions)"))
    end

    helper_D = copy(D)
    chosen_indices = Array{Tuple}(undef, number_of_reactions)
    for i in 1:number_of_reactions
        choice = rand(rng, non_zero_pos(helper_D))
        chosen_indices[i] = choice
        helper_D[choice[1], choice[2]] = 0
    end

    return chosen_indices
end

function checks_before_run(D, W_ba, tau, alpha)
    if isnothing(D)
        if !isnothing(W_ba)
            @warn "WARNING: Supplied energy yield matrix (W_ba) but no stoichiometric matrix (D). Creating D matrix of same size"
            D, W_ba = create_metabolism(n_resources=length(W_ba[1, :]))
        else
            D, W_ba = create_metabolism()
        end
    else
        if isnothing(W_ba)
            D, W_ba = create_metabolism(n_resources=length(D[1, :]))
            @warn "WARNING: Supplied stoichiometric matrix (D) but no energy yield matrix (W_ba). Creating W_ba matrix of same size"
        end
    end

    n_resources = size(D, 1)
    D_row, D_col = size(D)
    
    if (n_resources != D_row) || (n_resources != D_col)
        throw(DomainError("Number of resources does not match the size of the stoichiometric matrix (D) \n 
        number of resources is $n_resources, sizes of D are ($D_row, $D_col)"))
    end

    W_row, W_col = size(W_ba)
    if (n_resources != W_row) || (n_resources != W_col)
        throw(DomainError("Number of resources does not match the size of the energy yield matrix (W_ba) \n 
        number of resources is $n_resources, sizes of D are ($W_row, $W_col)"))
    end

    if isnothing(tau)
        tau = ones(Float64, n_resources)
    end

    if isnothing(alpha)
        alpha = vcat([100.0], zeros(Float64, n_resources-1))
    end

    if length(tau) != n_resources || length(alpha) != n_resources
        throw(DomainError("Length of dilution terms (tau) or resource availabilities (alpha) does not match the number of resources \n This can happen when tau or alpha is supplied, but the stoichiometric matrix and energy yield matrix are not. \n"))
    end

    return n_resources, D, W_ba, tau, alpha
end

function bounded_rand(rng, dist, bounds)
    out = rand(rng, dist)
    if out > bounds[1] && out < bounds[2]
        return out
    else
        return bounded_rand(rng, dist, bounds)
    end
end

function get_byproducts(monomer_content, byproduct_part, rng, allowed)
    byproducts = zeros(length(monomer_content))
    extra_splits = -1
    while byproduct_part >= 1
        extra_splits += 1
        possible_byproducts = findall(x -> x <= byproduct_part, monomer_content)
        for idx in possible_byproducts
            if allowed[idx] != 1
                possible_byproducts = filter(x -> x != idx, possible_byproducts)
            end
        end
        idx = rand(rng, possible_byproducts)
        byproduct_part -= monomer_content[idx]
        byproducts[idx] += 1
    end
    return byproducts, extra_splits
end

function partition_resources(n_resources, consumed_indices, consumption_rate, monomer_content, limiting_matrix, avg_usage, sd_usage, avg_efficiency, resource_sd, rng)
    consumption_rates = zeros(Float64, n_resources)
    energy = zeros(Float64, n_resources)
    production = zeros(Float64, n_resources, n_resources)
    n_splits = 0

    for idx in consumed_indices
        consumption_rates[idx] = consumption_rate
        energy_part = monomer_content[idx] * bounded_rand(rng, Normal(avg_usage, sd_usage), (0, 1))
        energy[idx] = energy_part * bounded_rand(rng, Normal(avg_efficiency, resource_sd[idx]), (0, 1))
        byproduct_part = floor(monomer_content[idx] - energy_part)
        if byproduct_part >= 1
            production[:, idx], extra_splits = get_byproducts(monomer_content, byproduct_part, rng, limiting_matrix[:, idx])
            n_splits += extra_splits
        end
        n_splits += energy_part
    end

    return consumption_rates, energy, production, n_splits
end

function shannon_diversity(abundances)
    total = sum(abundances)
    if total == 0 || isempty(abundances)
        return 0.0
    end
    proportions = abundances ./ total
    H = -sum(p * log(p) for p in proportions if p > 0)
    return H
end

function shannons_without_typing(data)
    shannon_vector = zeros(size(data)[2])
    for sample in 1:size(data)[2]
        shannon_vector[sample] = shannon_diversity(data[:, sample])
    end
    return shannon_vector
end

function describe_metab_matrix(D)
    n_consumed = dropdims(sum(!iszero, D, dims=1), dims=1)
    n_reactions = sum(!iszero, D)
    return n_consumed, n_reactions
end

function fuse_metabolism(carbon1, limiting1, carbon2, limiting2)
    # carbon1 and limiting1 define mucus metabolism
    size1 = length(carbon1)
    size2 = length(carbon2)
    orig_mucus = carbon1[1]
    carbon1[1] = maximum(carbon2) + 1

    full_limiting = zeros(size1 + size2, size1 + size2)
    full_limiting[1:size1, 1:size1] = limiting1
    full_limiting[size1+1:end, size1+1:end] = limiting2

    perm = sortperm(vcat(carbon1, carbon2), rev=true)
    full_limiting = full_limiting[:, perm]
    full_limiting = full_limiting[perm, :]
    full_carbon = sort(vcat(carbon1, carbon2), rev=true)
    full_carbon[1] = orig_mucus

    return full_carbon, Int.(full_limiting)
end