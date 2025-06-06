function describe_full_comm_species(full_ts)
    sizers = size(assay(full_ts[1], "non_negative"))
    total_abundances = zeros(Float64, sizers)

    alphadivs = [[] for i in keys(full_ts)]

    comm_tot_abundances_ts = []
    gen_tot_ts = []
    spec_tot_ts = []
    for i in keys(full_ts)
        alphadivs[i] = shannon(full_ts[i], "non_negative")
        total_abundances .+= assay(full_ts[i], "non_negative")
        comm_tot_ts = sum(assay(full_ts[i], "non_negative"), dims=1)
        push!(comm_tot_abundances_ts, vec(comm_tot_ts))
        gs_ts = aggregate_by_rowdata(full_ts[i], "n_reactions", "non_negative")
        gen_spec_keys = collect(keys(gs_ts))
        if length(gen_spec_keys) != 2
            raise("Error: Expected 2 keys in gen_spec_keys, got $(length(gen_spec_keys))")
        end
        gen_key = maximum(gen_spec_keys)
        spec_key = minimum(gen_spec_keys)

        if isempty(gen_tot_ts)
            append!(gen_tot_ts, gs_ts[gen_key])
            append!(spec_tot_ts, gs_ts[spec_key])
        else
            gen_tot_ts .+= gs_ts[gen_key]
            spec_tot_ts .+= gs_ts[spec_key]
        end
    end

    species_counts_ts = sum(!iszero, total_abundances, dims=1)
    println("length spec_counts: ", length(species_counts_ts))
    present_species_ts = [findall(!iszero, total_abundances[:, i]) for i in eachindex(total_abundances[1, :])]
    gammadivs = shannons_without_typing(total_abundances)

    betadivs = []
    for i in eachindex(total_abundances[1, :])
       alpha_mean = mean([alphadivs[j][i] for j in keys(full_ts)])
       append!(betadivs, gammadivs[i] / alpha_mean - 1)
    end

    return comm_tot_abundances_ts, vec(species_counts_ts), present_species_ts,
           gammadivs, betadivs, gen_tot_ts, spec_tot_ts
end


"""
This function takes an array of length equal to the number of communities.
For each community, an array of length equal to the number of timepoints is provided.
Each element of this array is an array of indices of species present at that timepoint in that community.

The function returns an array of arrays which describes the total species pool at each timepoint accross all communities.
"""
function total_species_pool(present_species_indices_ts::Array{Array{Array{Int64}, 1}, 1})
    ncomms = length(present_species_indices_ts)
    n_timepoints = length(present_species_indices_ts[1])
    total_pool = zeros(Int64, n_timepoints)
    for i in 1:n_timepoints
        total_pool[i] = collect(Set(Iterators.flatten([present_species_indices_ts[j][i] for j in 1:ncomms])))
    end
end

"""
This function returns the production matrices of the species present at each timepoint.
"""
function get_present_production_matrices(present_species_ts, all_matrices)
    present_matrices = []
    for present_species in present_species_ts
        matrices = all_matrices[present_species]
        push!(present_matrices, matrices)
    end
    return present_matrices
end


"""
This function takes the present production matrices at each timepoint
and calculates the number of different reactions (matrix positions) that are not zero.
"""
function total_metabolism(present_production_matrices_ts)
    n_timepoints = length(present_production_matrices_ts)
    println("n_timepoints:", n_timepoints)
    metabolic_capacity_ts = []
    for i in 1:n_timepoints
        total_metabolism = zeros(size(present_production_matrices_ts[i][1]))
        for matrix in present_production_matrices_ts[i]
            total_metabolism .+= matrix
        end
        append!(metabolic_capacity_ts, sum(total_metabolism .!= 0))
    end
    println("size mci:", size(metabolic_capacity_ts))
    return metabolic_capacity_ts
end