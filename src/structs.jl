
mutable struct ParamStruct
    n_species::Int64
    n_resources::Int64
    present_species::Array{Int64}
    consumption_rates::Array{Array{Float64}}
    production_matrices::Array{Union{Matrix{Int64}, Matrix{Float64}}}
    energy::Array{Array{Float64}}
    types::Array{String}
    n_reactions::Array{Int64}
    n_splits::Array{Float64}
    m::Array{Float64}
    phi::Float64
    eta::Float64
    tau::Array{Float64}
    alpha::Array{Float64}
    D_species::Float64
    D_resources::Float64
    A_species::Float64
    A_resources::Float64
    ph_opts::Array{Float64}
    ph_strength::Float64
    R_0::Float64
    decomp_percentage::Float64
    consumption_fluxes
    production_fluxes

    function ParamStruct(
        n_species::Int64, n_resources::Int64, 
         present_species::Array{Int64}, consumption_rates::Array{Array{Float64}}, production_matrices::Array{Union{Matrix{Int64}, Matrix{Float64}}}, 
          energy::Array{Array{Float64}}, types::Array{String}, n_reactions::Array{Int64}, n_splits::Array{Float64}, m::Array{Float64}, phi::Float64, eta::Float64, tau::Array{Float64}, 
            alpha::Array{Float64}; D_species::Float64=0.0, D_resources::Float64=0.0, A_species::Float64=0.0,
             A_resources::Float64=0.0, ph_opts::Array{Float64}=fill(7.0, n_species), ph_strength::Float64=1.0, R_0::Float64=100.0, decomp_percentage::Float64=0.0,
              consumption_fluxes=[], production_fluxes=[])
           
        new(n_species, n_resources, present_species, consumption_rates, production_matrices, energy, types, n_reactions, n_splits, m, phi, eta, tau, alpha, D_species, D_resources,
            A_species, A_resources, ph_opts, ph_strength, R_0, decomp_percentage, consumption_fluxes, production_fluxes)
    end

end

struct PoolStruct
    consumption_rates::Array{Array{Float64}}
    energy::Array{Array{Float64}}
    production_matrices::Array{Union{Matrix{Int64}, Matrix{Float64}}}
    m::Array{Float64}
    n_reactions::Array{Int64}
    n_splits::Array{Float64}
    types::Array{String}
end

struct SampleStruct
    n_species::Int64
    n_invaders::Int64
    consumption_rates::Array{Array{Float64}}
    energy::Array{Array{Float64}}
    production_matrices::Array{Union{Matrix{Int64}, Matrix{Float64}}}
    m::Array{Float64}
    n_reactions::Array{Int64}
    n_splits::Array{Float64}
    species_abundance::Array{Float64}
    resource_abundance::Array{Float64}
    ph_opts::Array{Float64}
    types::Array{String}
end