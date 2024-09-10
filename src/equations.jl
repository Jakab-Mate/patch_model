function equations(u, p, t; ph=7.0)
    species = zeros(Float64, p.n_species)
    consumption = zeros(Float64, p.n_resources)
    production = zeros(Float64, p.n_resources)

    for i in p.present_species
        species[i] = (sum(p.C[:,:,i] .* p.W_ba .* u[p.n_species+1:end]') * exp(-1 * (ph - p.ph_opts[i])^2) - (p.m[i] + p.phi * p.n_reactions[i] + p.eta * p.n_splits[i])) * u[i]
        if p.host_regulation
            species[i] = species[i] * (1 / (1 + exp(p.a[i] * (u[i]-p.k[i]))))
        end
        consumption += dropdims(sum(p.C[:, :, i] .* u[i] .* u[p.n_species+1:end]', dims=1)', dims=2)
        production += dropdims(sum(p.D .* p.C[:, :, i] .* u[i] .* u[p.n_species+1:end]', dims=2), dims=2)
    end

    resources = (p.alpha .- u[p.n_species+1:end]) ./ p.tau - consumption + production
    return vcat(species, resources)
end