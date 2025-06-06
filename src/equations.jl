"""
Calculate derivatives for each species and resource in the system.
"""
function equations(u, p, t; ph=7.0, alpha=nothing)
    species = zeros(Float64, p.n_species)
    consumption = zeros(Float64, p.n_resources)
    production = zeros(Float64, p.n_resources)
    for i in p.present_species
        species[i] = (sum(p.energy[i] .* u[p.n_species+1:end]) * exp(-p.ph_strength * (ph - p.ph_opts[i])^2) - (p.m[i] + p.eta * p.n_reactions[i] + p.phi * p.n_splits[i])) * u[i]
        consumption += (p.consumption_rates[i] .* u[i] .* u[p.n_species+1:end])
        production += dropdims(sum(p.production_matrices[i] .* u[i] .* u[p.n_species+1:end]', dims=2), dims=2)
    end
    if isnothing(alpha)
        resources = (p.alpha .- u[p.n_species+1:end-1]) ./ p.tau - consumption + production
    else
        resources = (alpha .- u[p.n_species+1:end]) ./ p.tau - consumption + production
    end
    return vcat(species, resources)
end


"""
Calculate derivatives for each species and resource in the system.
This version is designed to be thread-safe, allowing for parallel computation.
"""
function threaded_equations(u, p, t; ph=7.0, alpha=nothing)
    species = zeros(Float64, p.n_species)
    consumption = zeros(Float64, p.n_resources)
    production = zeros(Float64, p.n_resources)

    # Thread-safe accumulators
    consumption_local = Array{Float64}(undef, nthreads(), p.n_resources)
    production_local = Array{Float64}(undef, nthreads(), p.n_resources)

    # Initialize thread-local arrays to zero
    @threads for tid in 1:nthreads()
        consumption_local[tid, :] .= 0
        production_local[tid, :] .= 0
    end

    # Compute contributions in parallel
    @threads for i in p.present_species
        species[i] = (sum(p.energy[i] .* u[p.n_species+1:end]) * exp(-p.ph_strength * (ph - p.ph_opts[i])^2) -
                     (p.m[i] + p.eta * p.n_reactions[i] + p.phi * p.n_splits[i])) * u[i]

        tid = threadid()  # Current thread ID
        consumption_local[tid, :] .+= p.consumption_rates[i] .* u[i] .* u[p.n_species+1:end]
        production_local[tid, :] .+= dropdims(sum(p.production_matrices[i] .* u[i] .* u[p.n_species+1:end]', dims=2), dims=2)
    end

    # Aggregate results from all threads
    for tid in 1:nthreads()
        consumption .+= consumption_local[tid, :]
        production .+= production_local[tid, :]
    end

    # Resource dynamics
    if isnothing(alpha)
        resources = (p.alpha .- u[p.n_species+1:end-1]) ./ p.tau - consumption + production
    else
        resources = (alpha .- u[p.n_species+1:end]) ./ p.tau - consumption + production
    end

    return vcat(species, resources)
end


"""
Calculate derivatives for each species and resource in the system.
This version is used to save fluxes for each resource-species interaction.
"""
function saving_equations(u, p, t; ph=7.0, alpha=nothing)
    species = zeros(Float64, p.n_species)
    consumption = zeros(Float64, p.n_resources)
    production = zeros(Float64, p.n_resources)
    for i in p.present_species
        species[i] = (sum(p.energy[i] .* u[p.n_species+1:end]) * exp(-p.ph_strength * (ph - p.ph_opts[i])^2) - (p.m[i] + p.eta * p.n_reactions[i] + p.phi * p.n_splits[i])) * u[i]
        consumption += (p.consumption_rates[i] .* u[i] .* u[p.n_species+1:end])
        push!(p.consumption_fluxes[1], p.consumption_rates[i] .* u[i] .* u[p.n_species+1:end])
        production += dropdims(sum(p.production_matrices[i] .* u[i] .* u[p.n_species+1:end]', dims=2), dims=2)
        push!(p.production_fluxes[1], dropdims(sum(p.production_matrices[i] .* u[i] .* u[p.n_species+1:end]', dims=2), dims=2))
    end
    if isnothing(alpha)
        resources = (p.alpha .- u[p.n_species+1:end-1]) ./ p.tau - consumption + production
    else
        resources = (alpha .- u[p.n_species+1:end]) ./ p.tau - consumption + production
    end
    return vcat(species, resources)
end

