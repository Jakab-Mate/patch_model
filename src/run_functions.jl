include("structs.jl")
include("plot_functions.jl")

function generic_run(path, D, W_ba, sample::sample_struct; 
    t_span::Tuple{Int64, Int64}=(0, 1000), 
    n_resources::Int64=10, 
    n_species::Int64=10, 
    n_invaders::Int64=0,
    t_inv::Float64=25.0, 
    t_inv_0::Float64=100.0, 
    cutoff::Float64=0.0001, 
    phi::Float64=0.1, 
    eta::Float64=0.1, 
    tau::Union{Vector{Float64}, Nothing}=nothing, 
    alpha::Union{Vector{Float64}, Nothing}=nothing,
    plot::Bool=true,
    host_regulation::Bool=true)

    if isnothing(tau)
        tau = ones(Float64, n_resources)
    end

    if isnothing(alpha)
        alpha = vcat([100.0], zeros(Float64, n_resources-1))
    end

    u0 = vcat(sample.species_abundance, sample.resource_abundance)
    params = param_struct(n_species+n_invaders, n_resources, 1:n_species, sample.C, D, W_ba, sample.n_reactions, sample.n_splits, sample.m, phi, eta, tau, alpha, sample.a, sample.k, host_regulation)

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
        condition = rand(["lake", "ocean", "river"], length(t)),
        time = 1:length(t)
    );

    se = SummarizedExperiment(assays, rowdata, coldata);
    if plot
        plot_se(se, "sim", path)
        plot_MCI(se, "sim", path)
    end
    return se

end

function sample_same_pool(n_samples, seed, rng, D, W_ba, pool::pool_struct; 
    t_span::Tuple{Int64, Int64}=(0, 1000), 
    n_resources::Int64=10, 
    n_species::Int64=10, 
    n_invaders::Int64=0,
    n_families::Int64=10,
    family_size::Int64=20, 
    t_inv::Float64=25.0, 
    t_inv_0::Float64=100.0, 
    cutoff::Float64=0.0001, 
    phi::Float64=0.1, 
    eta::Float64=0.1, 
    tau::Union{Vector{Float64}, Nothing}=nothing, 
    alpha::Union{Vector{Float64}, Nothing}=nothing)

    if isnothing(tau)
        tau = ones(Float64, n_resources)
    end
    if isnothing(alpha)
        vcat([100.0], zeros(Float64, n_resources-1))
    end

    for s in 1:n_samples
        sample = sample_pool(rng, pool, n_resources, n_species, n_invaders)
        out = generic_run(path, D, W_ba, sample,
            t_span=t_span, 
            n_resources=n_resources, 
            n_species=n_species, 
            n_invaders=n_invaders, 
            t_inv=t_inv, 
            t_inv_0=t_inv_0, 
            cutoff=cutoff, 
            phi=phi, 
            eta=eta, 
            tau=tau,
            alpha=alpha)

        # add to ordered dict --- assays = OrderedDict{String, AbstractArray}("sim" => Xapp)
        
    end
end

