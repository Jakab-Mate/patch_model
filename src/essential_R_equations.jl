"""
Calculate derivatives for each species and resource in the system.
This version includes an essential resource and its effects on species.
"""
function essential_R_equations(u, p::ParamStruct, t; ph=7.0, alpha=nothing)
    species = zeros(Float64, p.n_species)
    consumption = zeros(Float64, p.n_resources)
    production = zeros(Float64, p.n_resources+1)
    
    E = u[end] / (p.R_0 + u[end])
    # E = 1 / sum(u[1:p.n_species]) * u[end] / (p.R_0 + u[end])
    total_decomp = 0

    for i in p.present_species
        total_decomp += (p.m[i] + p.eta * p.n_reactions[i] + p.phi * p.n_splits[i]) * u[i]
        species[i] = (E * sum(p.energy[i] .* u[p.n_species+1:end-1]') * exp(-p.ph_strength * (ph - p.ph_opts[i])^2) - (p.m[i] + p.eta * p.n_reactions[i] + p.phi * p.n_splits[i])) * u[i]
        consumption += dropdims((p.consumption_rates[i] .* u[i] .* u[p.n_species+1:end]')', dims=2)
        production += dropdims(sum(p.production_matrices[i] .* u[i] .* u[p.n_species+1:end]', dims=2), dims=2)
    end

    if isnothing(alpha)
        regular_resources = (p.alpha .- u[p.n_species+1:end-1]) ./ p.tau - consumption + production
    else
        regular_resources = (alpha .- u[p.n_species+1:end-1]) ./ p.tau - consumption + production
    end
    essential_resource = p.decomp_percentage * total_decomp - E * sum(u[1:p.n_species])
    return vcat(species, regular_resources, essential_resource)
end