"""
    generic_run(sample::SampleStruct;
       D=nothing,
       W_ba=nothing,
       path=homedir(),
       t_span=(0, 1000),
       t_inv=25.0,
       t_inv_0=100.0,
       cutoff=0.0001,
       phi=0.1,
       eta=0.1,
       tau=nothing,
       alpha=nothing,
       plot=true,
       host_regulation=true)

Run the model with the given parameters and sample by solving a set of ODEs using the KenCarp4 solver.

# Mandatory arguments
- `sample::SampleStruct`: A sample struct containing the initial conditions and parameters for the model.

# Recommended but optional arguments
- `D::Union{Nothing, Matrix{Int64}}`: The stoichiometric matrix for the model. If not supplied, a default matrix will be created.
- `W_ba::Union{Nothing, Matrix{Float64}}`: The energy yield matrix for the model. If not supplied, a default matrix will be created.
- `path::String`: The path to save the output plots. Default is `homedir()`.
- `t_span::Tuple{Int64, Int64}`: The time span for the simulation. Default is `(0, 1000)`.

# Optional arguments
- `t_inv::Float64`: The time between the introduction of two subsequent invading species. Default is `25.0`.
- `t_inv_0::Float64`: The time at which the first invading species is introduced. Default is `100.0`.
- `cutoff::Float64`: The abundance threshold under which a species is considered extinct and its abundance is set to 0. Default is `0.0001`.
- `phi::Float64`: The strength of the additional maintenance costs based on the complexity of a the reaction repertoires of species. Default is `0.1`.
- `eta::Float64`: The strength of the additional maintenance costs based on the number of reactions of species. Default is `0.1`.
- `tau::Union{Vector{Float64}, Nothing}`: Controls the replenisment/depletion rates of resources from/into the outter environment. Default is `1.0` for all reasources.
- `alpha::Union{Vector{Float64}, Nothing}`: The availability of resources in the outer environment. Default is `100.0` for the first resource and 0.0 for the rest.
- `plot::Bool`: Whether to generate plots of the simulation. Default is `true`.
- `host_regulation::Bool`: Whether to include host regulation in the model. Default is `true`.

# Output
Returns time series data in a SummarizedExperiment (SE) data container, which can be used for a variety of analyses. For details, see [MicrobiomeAnalysis.jl](https://github.com/JuliaTurkuDataScience/MicrobiomeAnalysis.jl)
"""
function generic_run(sample::SampleStruct;
    D::Union{Nothing, Matrix{Int64}, Matrix{Float64}}=nothing,
    W_ba::Union{Nothing, Matrix{Float64}}=nothing,
    path::String=homedir(),
    t_span::Tuple{Int64, Int64}=(0, 1000), 
    t_inv::Float64=25.0, 
    t_inv_0::Float64=100.0, 
    cutoff::Float64=0.0001, 
    phi::Float64=0.1, 
    eta::Float64=0.1, 
    tau::Union{Vector{Float64}, Nothing}=nothing, 
    alpha::Union{Vector{Float64}, Nothing}=nothing,
    plot::Bool=true,
    decomp_percentage::Float64=0.1)

    n_resources, D, W_ba, tau, alpha = checks_before_run(D, W_ba, tau, alpha)

    n_species = sample.n_species
    n_invaders = sample.n_invaders

    u0 = vcat(sample.species_abundance, sample.resource_abundance)

    consumption_rates = []
    production_matrices = []
    energy = []
    for C in eachslice(sample.C, dims=3)
        push!(consumption_rates, sum(C[1:size(W_ba, 1), 1:size(W_ba, 2)], dims=1))
        push!(production_matrices, D .* C)
        push!(energy, sum(W_ba .* C[1:size(W_ba, 1), 1:size(W_ba, 2)], dims=1))
    end

    println(typeof(energy))
    println(typeof(consumption_rates))
    println(typeof(production_matrices))

    params = ParamStruct(n_species+n_invaders, n_resources, collect(1:n_species), consumption_rates, production_matrices, energy, sample.n_reactions,
        sample.n_splits, sample.m, phi, eta, tau, alpha, sample.a, sample.k, decomp_percentage=decomp_percentage)



    prob = ODEProblem(equations, u0, t_span, params)

    
    cb = create_callbacks(t_inv, t_inv_0, n_invaders, n_species, cutoff)
    solution = solve(prob, KenCarp4(autodiff=false), abstol = 1e-10, reltol = 1e-10; saveat=1, callback=cb)

    assay_name = "sim"
    t = solution.t
    C_vec = [sample.C[:, :, i] for i in 1:size(sample.C, 3)]
    println("Size of solution.u: ", size(solution.u))
    Xapp = hcat(map(u -> u[1:(n_species+n_invaders)], solution.u)...)
    assays = OrderedDict{String, AbstractArray}(assay_name => Xapp);
    rowdata = DataFrame(
        name = ["strain$i" for i in 1:(n_species+n_invaders)],
        family = ["family$i" for i in sample.family_ids],
        n_reactions = ["$i reactions" for i in sample.n_reactions],
        n_splits = ["$i splits" for i in sample.n_splits],
        ph_opt = ["pH $i" for i in sample.ph_opts],
        matrix = C_vec
    );
    coldata = DataFrame(
        name = ["t$i" for i in 1:length(t)],
        time = 1:length(t)
    );

    se = SummarizedExperiment(assays, rowdata, coldata);
    if plot
        plot_se(se, assay_name, path)
        plot_MCI(se, assay_name, path)
    end
    return se, assay_name # Document! Now returning assay name as well.

end

"""
    spatial_run(n_comms::Int64, sample::SampleStruct;
        ph_list::Union{Nothing, Vector{Float64}}=nothing,
        ph_strength=1.0,
        D_species=0.1,
        D_resources=0.1,
        A_species=0.1,
        A_resources=0.1,
        D=nothing,
        W_ba=nothing,
        path=homedir(),
        t_span=(0, 1000),
        t_inv=25.0,
        t_inv_0=100.0,
        cutoff=0.0001,
        phi=0.1,
        eta=0.1,
        tau=nothing,
        alpha=nothing,
        plot=true,
        host_regulation=false)

Run the spatially extended version of the model, where a set of communities are simulated. The communities resemble different sections of the gut,
where resources and invaders appear in the first community and can spread to the other communities through diffusion and advection.

# Mandatory arguments
- `n_comms::Int64`: The number of communities to simulate.
- `sample::SampleStruct`: A sample struct containing the initial conditions and parameters for the model.

# Additional arguments
- `D_species::Float64`: The diffusion coefficient for species. Default is `0.1`.
- `D_resources::Float64`: The diffusion coefficient for resources. Default is `0.1`.
- `A_species::Float64`: The advection coefficient for species. Default is `0.1`.
- `A_resources::Float64`: The advection coefficient for resources. Default is `0.1`.
- `ph_strength::Float64`: The strength of the pH effect on growth. Default is `1.0`.
- `ph_list::Union{Nothing, Vector{Float64}}`: A list of pH values for each community. If not supplied, the default pH value is 7.0.

# Output
Returns a dictionary of SummarizedExperiment (SE) data containers
"""
function spatial_run(n_comms::Int64, sample::SampleStruct, seed::Int64;
    ph_list::Union{Nothing, Vector{Float64}} = nothing,
    ph_strength::Float64=1.0,
    D_species::Float64=0.1,
    D_resources::Float64=0.1,
    A_species::Float64=0.1,
    A_resources::Float64=0.1,
    path::String=homedir(),
    t_span::Tuple{Int64, Int64}=(0, 1000),
    t_inv::Float64=25.0,
    t_inv_0::Float64=1.0,
    cutoff::Float64=0.0001,
    phi::Float64=0.1,
    eta::Float64=0.1,
    tau::Union{Nothing, Array{Float64}}=nothing,
    alpha::Union{Nothing, Array{Float64}}=nothing,
    plot::Bool=true,
    alpha_list::Union{Nothing, Vector{Any}}=nothing,
    symmetrical::Bool=false,
    plotgraph::Bool=false,
    consume_all_monomers::Int64=0,
    influx_amount::Float64=0.0,
    influx_time::Float64=0.0,
    n_complex::Int64=3,
    mucus::Bool=false)

    rng = MersenneTwister(seed)

    if influx_amount > 0.0
        periodic_influx = true
    else
        periodic_influx = false
    end

    n_resources = length(sample.production_matrices[1][1, :])

    if isnothing(alpha)
        alpha = zeros(Float64, n_resources)
        alpha[1] = 100.0
    end

    if isnothing(tau)
        tau = ones(Float64, n_resources)
    end

    if isnothing(alpha_list)
        alpha_list = []
        push!(alpha_list, alpha)
        for i in 1:n_comms-1
            push!(alpha_list, zeros(Float64, n_resources))
        end
    end

    n_species = sample.n_species
    n_invaders = sample.n_invaders
    ####

    u0 = vcat(sample.species_abundance, sample.resource_abundance)

    params = ParamStruct(n_species+n_invaders, n_resources, collect(1:n_species), sample.consumption_rates, sample.production_matrices, sample.energy, 
            sample.types, sample.n_reactions, sample.n_splits, sample.m, phi, eta, tau, alpha, D_species=D_species, D_resources=D_resources, 
            A_species=A_species, A_resources=A_resources, ph_opts=sample.ph_opts, ph_strength=ph_strength)

    u_mat = zeros(Float64, n_species+n_invaders+n_resources, n_comms)
    du_mat = zeros(Float64, n_species+n_invaders+n_resources, n_comms)


    du = zeros(Float64, length(u0))
    
    prob = ODEProblem((du, u, p, t) -> spatial_equations!(du, u, p, t, ph_list=ph_list, N=n_comms, u_mat=u_mat, du_mat=du_mat, alpha_list=alpha_list, symmetrical=symmetrical),
            u0, t_span, params)

    invasion_cb, resource_cb = create_spatial_callbacks(t_inv, t_inv_0, n_invaders, n_species, cutoff, n_comms, influx_amount, influx_time, n_complex, mucus, rng)
    if periodic_influx
        cb = CallbackSet(invasion_cb, resource_cb)
    else
        cb = CallbackSet(invasion_cb)
    end
    println("Starting simulation")
    solution = solve(prob, KenCarp4(autodiff=false), abstol = 1e-10, reltol = 1e-10; saveat=t_inv/50, callback=cb)
    println("Simulation finished")

    if plotgraph
        solution.prob.p.consumption_fluxes = Vector{Vector{Vector{Float64}}}(undef, n_comms + 1)
        println("Size of consumption fluxes: ", size(solution.prob.p.consumption_fluxes))
        solution.prob.p.production_fluxes = Vector{Vector{Vector{Float64}}}(undef, n_comms + 1)
        println("Size of production fluxes: ", size(solution.prob.p.production_fluxes))

        restart_prob = ODEProblem((du, u, p, t) -> saving_spatial_equations!(du, u, p, t, N=n_comms, u_mat=u_mat, du_mat=du_mat, alpha_list=alpha_list, symmetrical=symmetrical),
                solution.u[end], (solution.t[end], t_span[2]+1), solution.prob.p)

        println("Restarting simulation")
        restart_sol = solve(restart_prob, KenCarp4(autodiff=false), abstol = 1e-10, reltol = 1e-10; saveat=1, callback=cb)

        species_abundances = []
        for i in 1:n_comms
            commspecies = solution.u[end][(i-1)*(n_species+n_invaders)+1:i*(n_species+n_invaders)]
            push!(species_abundances, commspecies[solution.prob.p.present_species])
        end

        plot_community_graph(sample.consumption_rates[solution.prob.p.present_species],
                             sample.production_matrices[solution.prob.p.present_species],
                             restart_sol.prob.p.consumption_fluxes[2:end],
                             restart_sol.prob.p.production_fluxes[2:end],
                             species_abundances,
                             path,
                             consume_all_monomers=consume_all_monomers)
    end

    t= solution.t
    se_dict_species = Dict{Int64, SummarizedExperiment}()
    se_dict_resources = Dict{Int64, SummarizedExperiment}()
    for i in 1:n_comms
        species_ts = hcat(map(u -> u[(i-1)*(n_species+n_invaders)+1:i*(n_species+n_invaders)], solution.u)...)
        non_negative = deepcopy(species_ts)
        for i in eachindex(non_negative)
            if non_negative[i] < cutoff
                non_negative[i] = 0.0
            end
        end
        species_assays = OrderedDict{String, AbstractArray}("species" => species_ts, "non_negative" => non_negative);
        species_rowdata = DataFrame(
            name = ["strain$i" for i in 1:(n_species+n_invaders)],
            n_reactions = ["$i reactions" for i in sample.n_reactions],
            n_splits = ["$i splits" for i in sample.n_splits],
            matrices = sample.production_matrices
        );
        coldata = DataFrame(
            name = ["t$i" for i in 1:length(t)],
            time = t
        );

        se_species = SummarizedExperiment(species_assays, species_rowdata, coldata);
        se_dict_species[i] = se_species

        r_idx0 = n_comms*(n_species+n_invaders) + 1
        resources_ts = hcat(map(u -> u[r_idx0+n_resources*(i-1):r_idx0+n_resources*i-1], solution.u)...)
        resource_assays = OrderedDict{String, AbstractArray}("resources" => resources_ts);
        resources_rowdata = DataFrame(
            name = ["resource$i" for i in 1:n_resources],
            number = [i for i in 1:n_resources]
        );

        se_resources = SummarizedExperiment(resource_assays, resources_rowdata, coldata)
        se_dict_resources[i] = se_resources
    end

    return se_dict_species, se_dict_resources 
end


# data_to_save = DataFrame(
#             name = ["strain$i" for i in 1:(n_species+n_invaders)],
#             n_reactions = ["$i reactions" for i in sample.n_reactions],
#             n_splits = ["$i splits" for i in sample.n_splits],
#         )

# for t_idx in 1:length(t)
#     data_to_save[!, Symbol("$t_idx")] = Xapp[:, t_idx]
# end

# #CSV.write(joinpath(path, "comm$i" * "_data.csv"), data_to_save)

# if plot
#     plot_se(se, "sim", path, comm=i)
#     plot_MCI(MCI_ts, n_species_ts, t)
# end