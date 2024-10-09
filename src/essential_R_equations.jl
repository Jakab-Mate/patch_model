function essential_R_equations(u, p, t; ph=7.0)
    species = zeros(Float64, p.n_species)
    consumption = zeros(Float64, p.n_resources)
    production = zeros(Float64, p.n_resources)

    E = 1 / sum(u[1:p.n_species]) * u[end] / (R_0 + u[end])
    total_decomp = 0

    for i in p.present_species
        total_decomp += (p.m[i] + p.eta * p.n_reactions[i] + p.phi * p.n_splits[i]) * u[i]
        species[i] = (sum(p.C[1:end-1, 1:end-1, i] .* p.W_ba .* u[p.n_species+1:end-1]') * exp(-p.ph_strength * (ph - p.ph_opts[i])^2) - (p.m[i] + p.eta * p.n_reactions[i] + p.phi * p.n_splits[i])) * u[i]
        if p.host_regulation
            species[i] = species[i] * (1 / (1 + exp(p.a[i] * (u[i]-p.k[i]))))
        end
        consumption += dropdims(sum(p.C[1:end-1, 1:end-1, i] .* u[i] .* u[p.n_species+1:end-1]', dims=1)', dims=2)
        production += dropdims(sum(p.D .* p.C[:, :, i] .* u[i] .* u[p.n_species+1:end]', dims=2), dims=2)
    end

    regular_resources = (p.alpha .- u[p.n_species+1:end]) ./ p.tau - consumption + production[1:end-1]
    essential_resource = decomp_percentage * total_decomp + production[end] - E * sum(u[1:p.n_species])
    resources = vcat(regular_resources, essential_resource)
    return vcat(species, resources)
end