### Functions that are not used

function describe_single_comm_species(ts, consumption_rates, production_matrices)
    n_consumed = []
    n_reactions = []
    for matrix in production_matrices
        cons, react = describe_metab_matrix(matrix)
        push!(n_consumed, cons)
        push!(n_reactions, react)
    end

    # Count species present at each timepoint
    species_counts_ts = sum(!iszero, ts, dims=2)

    # Find indices of species present at each timepoint
    species_indices_ts = [findall(!iszero, col) for col in eachcol(A)]

    # Find production matrices of species present at each timepoint
    production_matrices_ts = [production_matrices[present] for present in species_indices_ts]

    # Find consumption rates of species present at each timepoint
    consumption_rates_ts = [consumption_rates[present] for present in species_indices_ts]

    # Find number of resources consumed by each species present at each timepoint
    n_consumed_ts = [n_consumed[present] for present in species_indices_ts]

    # Find number of reactions of each species present at each timepoint
    n_reactions_ts = [n_reactions[present] for present in species_indices_ts]

    return species_counts_ts, species_indices_ts,
           production_matrices_ts, consumption_rates_ts,
           n_consumed_ts, n_reactions_ts
end

"""
Remove any species with abundance below `cutoff` from the list of present species.
If `start_time` has passed and the time since the last invader addition is greater than `t_inv`,
add an invader, by setting its previously 0 abundance to 10.
"""
function affect!(integrator, t_inv, start_time, n_invaders, n_species, cutoff)
    integrator.p.present_species = findall(x -> x > cutoff, integrator.u[1:(n_species+n_invaders)])
    if integrator.t >= start_time
        idx = Int(n_species + floor((integrator.t - start_time) / t_inv)) + 1
        if idx <= n_invaders + n_species
            println("invader ", idx-n_species, " added at time ", integrator.t)
            integrator.u[idx] += 10
        end
    end
end