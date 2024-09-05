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
- `D::Union{Nothing, AbstractMatrix}`: The stoichiometric matrix for the model. If not supplied, a default matrix will be created.
- `W_ba::Union{Nothing, AbstractMatrix}`: The energy yield matrix for the model. If not supplied, a default matrix will be created.
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
    D = nothing,
    W_ba = nothing,
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
    host_regulation::Bool=true)

    D, W_ba, tau, alpha = checks_before_run(D, W_ba, tau, alpha)

    n_species = sample.n_species
    n_invaders = sample.n_invaders

    u0 = vcat(sample.species_abundance, sample.resource_abundance)
    params = ParamStruct(n_species+n_invaders, n_resources, collect(1:n_species), sample.C, D, W_ba, sample.n_reactions, sample.n_splits, sample.m, phi, eta, tau, alpha, sample.a, sample.k, host_regulation)

    prob = ODEProblem(equations, u0, t_span, params)
    cb = create_callbacks(t_inv, t_inv_0, n_invaders, n_species, cutoff)
    solution = solve(prob, KenCarp4(autodiff=false), abstol = 1e-10, reltol = 1e-10; saveat=1, callback=cb)

    t = solution.t
    Xapp = hcat(map(u -> u[1:(n_species+n_invaders)], solution.u)...)
    assays = OrderedDict{String, AbstractArray}("sim" => Xapp);
    rowdata = DataFrame(
        name = ["strain$i" for i in 1:(n_species+n_invaders)],
        family = ["family$i" for i in sample.family_ids],
        n_reactions = ["$i reactions" for i in sample.n_reactions],
        n_splits = ["$i splits" for i in sample.n_splits]
    );
    coldata = DataFrame(
        name = ["t$i" for i in 1:length(t)],
        time = 1:length(t)
    );

    se = SummarizedExperiment(assays, rowdata, coldata);
    if plot
        plot_se(se, "sim", path)
        plot_MCI(se, "sim", path)
    end
    return se

end

"""
    spatial_run(n_comms::Int64, sample::SampleStruct;
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
    host_regulation=true)

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

# Output
Returns a dictionary of SummarizedExperiment (SE) data containers
"""
function spatial_run(n_comms::Int64, sample::SampleStruct;
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
    host_regulation=true)

    n_resources, D, W_ba, tau, alpha = checks_before_run(D, W_ba, tau, alpha)

    n_species = sample.n_species
    n_invaders = sample.n_invaders

    u0 = vcat(vcat(sample.species_abundance, repeat(zeros(n_species+n_invaders), n_comms-1)), vcat(sample.resource_abundance, repeat(zeros(n_resources), n_comms-1)))

    params = ParamStruct(n_species+n_invaders, n_resources, collect(1:n_species), sample.C, D, W_ba, sample.n_reactions, sample.n_splits, sample.m,
            phi, eta, tau, alpha, sample.a, sample.k, host_regulation, D_species=D_species, D_resources=D_resources, A_species=A_species, A_resources=A_resources)

    prob = ODEProblem(spatial_equations, u0, t_span, params)
    cb = create_spatial_callbacks(t_inv, t_inv_0, n_invaders, n_species, cutoff, n_comms)
    solution = solve(prob, KenCarp4(autodiff=false), abstol = 1e-10, reltol = 1e-10; saveat=1, callback=cb)

    t= solution.t
    se_dict = Dict{Int64, SummarizedExperiment}()
    for i in 1:n_comms
        Xapp = hcat(map(u -> u[(i-1)*(n_species+n_invaders)+1:i*(n_species+n_invaders)], solution.u)...)
        assays = OrderedDict{String, AbstractArray}("sim" => Xapp);
        rowdata = DataFrame(
            name = ["strain$i" for i in 1:(n_species+n_invaders)],
            family = ["family$i" for i in sample.family_ids],
            n_reactions = ["$i reactions" for i in sample.n_reactions],
            n_splits = ["$i splits" for i in sample.n_splits]
        );
        coldata = DataFrame(
            name = ["t$i" for i in 1:length(t)],
            time = 1:length(t)
        );

        se = SummarizedExperiment(assays, rowdata, coldata);
        se_dict[i] = se
        if plot
            plot_se(se, "sim", path, comm=i)
            plot_MCI(se, "sim", path, comm=i)
        end
    end

    return se_dict

end