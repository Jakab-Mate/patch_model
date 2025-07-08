"""
MANDATORY ARGUMENTS
- n_repeats::Int64: Number of repetitions for the simulation.
- basepath::String: Path to the directory where the results will be saved.

OPTIONAL ARGUMENTS  
Resource parameters:
- `n_complex`::Int64: Number of complex resources. Default is `3`.
- `n_simple`::Int64: Number of simple resources. Default is `10`.
- `n_monomer`::Int64: Number of monomer resources. Default is `5`.
- `gapsize`::Int64: Size of the gap between complex and simple resources. Default is `15`.

Resource supply parameters:
- `influx_amount`::Float64: Amount of resource influx. Default is `0.0`.
- `influx_time`::Float64: Time between resource influxes. Default is `1.0`.
- `alpha_high`::Float64: Controls constant resource supply if `mucus=false`. Each local community receives a continuous supply 100 units of complex resource, of which `alpha_high` amount is dedicated to a single complex resource. The remainder is split evenly between other complex resources. Default is `100.0`.

Species parameters:
- `n_invaders`::Int64: Number of invading species. Default is `100`.
- `t_inv`::Float64: Time between the introduction of invading species. Default is `1.0`.
- `n_species`::Int64: Number of species initially in the habitat. Default is `1`.
- `pool_size`::Int64: Size of the species pool. Default is `500`.
- `tradeoff_strength`__number_consumed::Float64: Trade-off strength for the number of resources consumed. Default is `0.1`.
- `tradeoff_strength`__number_of_splits::Float64: Trade-off strength for the "complexity" of a species' reaction repertoire, based on the number of enzymatic cuts they perform. Default is `0.1`.
- `h`::Float64: Reaction rate allocation tradeoff. The h-th power of reaction rates in a species always sums to one. Default is `1.0`.
- `gen_n_consumed`::Int64: Number of resources consumed by generalist species. Default is `4`.
- `spec_n_consumed`::Int64: Number of resources consumed by specialist species. Default is `1`.
- `consume_all_monomers`::Bool: Whether every species is able to consume all monomers. Default is `false`.

Mucus parameters:
- `mucus`::Bool: Whether to include mucus in the simulation, where resource supply is constant. Default is `true`.
- `mucus_disjunct`::Bool: Whether the mucus resource is disjunct from the rest of the primary resources (there is no overlap in their byproducts). Default is `true`.
- `mucus_simple`::Int64: Only used if `mucus_disjunt=true` Number of simple resources in the mucus resource breakdown pathways. Default is `5`.
- `mucus_mono`::Int64: Only used if `mucus_disjunt=true` Number of monomer resources in the mucus resource breakdown pathways. Default is `4`.
- `mucus_gap`::Int64: Only used if `mucus_disjunt=true` Size of the gap between mucus and simple mucus-specific resources. Default is `10`.

Simulation parameters:
- `ncomms`::Int64: Number of communities in the simulation. Default is `3`.
- `t_span`::Tuple{Float64, Float64}: Time span of the simulation. Default is `(0, 150)`.
- `cutoff`::Float64: Species under cutoff are considered extinct. Default is `0.0001`.
- `symmetrical`::Bool: Whether diffusion is symmetrical or directional across communities. Default is `true`.
- `limited_pathways`::Bool: Whether to limit the underlying metabolic reaction system by excluding ~half of the energetically possible byproducts of resources. Default is `true`.

Stochasticity parameters:
- `set_seeds`::Union{Nothing, Vector{Int64}}: If provided, seeds for each repeat. Default is `nothing`.
- `fix_metab`::Bool: Whether to use the same metabolism for all repeats. Default is `false`.
- `fix_pool`::Bool: Whether to use the same species pool for all repeats. Default is `false`.
- `fix_sample`::Bool: Whether to use the same sample for all repeats. Default is `false`.
- `fix_callbacks`::Bool: Whether to use the same callbacks for all repeats. Default is `false`.
- `fix_order`::Bool: Whether to fix the order of invading species. Default is `false`.
"""
function repeat_params(n_repeats, basepath; 
    # Resource params
    n_complex=3,
    n_simple=10,
    n_monomer=5,
    gapsize=15,
    n_resources = n_complex + n_simple + n_monomer,
    # Resource supply params
    influx_amount = 0.0,
    influx_time = 1.0,
    alpha_high = 100.0, #TODO this should be a vector for more customization
    # Species params
    n_species = 1,
    n_invaders = 100,
    t_inv = 1.0,
    pool_size = 500,
    tradeoff_strength__number_consumed = 0.1,
    tradeoff_strength__number_of_splits = 0.1,
    h = 1.0,
    gen_n_consumed = 4,
    spec_n_consumed = 1,
    consume_all_monomers = true,

    # Mucus params
    mucus = true,
    mucus_disjunct = true,
    mucus_simple = 5,
    mucus_mono = 4,
    mucus_gap = 10,
    # Simulation params
    ncomms = 3,
    t_span = (0, 150),
    cutoff = 0.0001,
    symmetrical = true,
    limited_pathways = true,
    set_seeds = nothing,
    fix_metab = false,
    fix_pool = false,
    fix_sample = false,
    fix_callbacks = false,
    fix_order = false)

    alpha_list = []
    if !mucus
        alpha_low = (100 - alpha_high) / n_complex
        for i in 1:ncomms
            alpha = vcat(fill(alpha_low, n_complex), zeros(Float64, n_resources-n_complex))
            if i == ncomms
                alpha[i] = alpha_high
            end
            push!(alpha_list, alpha)
        end
    else
        n_resources = n_resources + mucus_simple + mucus_mono
        alpha = zeros(Float64, n_resources)
        alpha[1] = alpha_high
        push!(alpha_list, alpha)
        for i in 2:ncomms
            push!(alpha_list, zeros(Float64, n_resources))
        end
    end

    shannons = [[], [], []]
    species_counts = [[], [], []]
    total_comm_tots = [[], [], []]
    total_shannons = []
    total_counts = []
    total_beta = []
    total_MCI = []
    total_gen = []
    total_spec = []
    seeds = []

    gen_tr_rel = [[], []]
    spec_tr_rel = [[], []]
    counts_tr_rel = [[], []]
    shannons_tr_rel = [[], []]
    MCI_tr_rel = [[], []]
    masterseed = trunc(Int, rand() * 10000)

    for i in 1:n_repeats
        if !isnothing(set_seeds)
            seed = set_seeds[i]
        else
            seed = trunc(Int, rand() * 10000)
        end
        append!(seeds, seed)
        path="$basepath/seed_$seed"
        mkpath(path)
        run_params = Dict("n_complex" => n_complex, "n_simple" => n_simple, "n_monomer" => n_monomer, "gapsize" => gapsize, "n_resources" => n_resources, "n_species" => n_species, "n_invaders" => n_invaders,
            "pool_size" => pool_size, "ncomms" => ncomms, "t_span" => t_span, "cutoff" => cutoff, "tradeoff_strength__number_consumed" => tradeoff_strength__number_consumed,
            "tradeoff_strength__number_of_splits" => tradeoff_strength__number_of_splits, "gen_n_consumed" => gen_n_consumed, "spec_n_consumed" => spec_n_consumed,
            "symmetrical" => symmetrical, "limited_pathways" => limited_pathways, "seed" => seed, "path" => path, "h"=>h,
            "consume_all_monomers" => consume_all_monomers, "influx_time" => influx_time, "influx_amount" => influx_amount, "mucus" => mucus, "t_inv" => t_inv,
            "mucus_disjunct" => mucus_disjunct, "mucus_simple" => mucus_simple, "mucus_mono" => mucus_mono, "mucus_gap" => mucus_gap)
        open("$path/params.txt", "w") do file
            for (key, value) in run_params
                println(file, "$key: $value")
            end
        end
        ########## SIMULATION ###########
        if mucus && mucus_disjunct
            if fix_metab
                m_carbon, m_limiting = create_metabolism(1, mucus_simple, mucus_mono, mucus_gap, limited_pathways=limited_pathways, seed=masterseed)
                reg_carbon, reg_limiting = create_metabolism(n_complex-1, n_simple, n_monomer, gapsize, limited_pathways=limited_pathways, seed=masterseed)
            else
                m_carbon, m_limiting = create_metabolism(1, mucus_simple, mucus_mono, mucus_gap, limited_pathways=limited_pathways, seed=seed)
                reg_carbon, reg_limiting = create_metabolism(n_complex-1, n_simple, n_monomer, gapsize, limited_pathways=limited_pathways, seed=seed)
            end
            carbon, limiting_matrix = fuse_metabolism(m_carbon, m_limiting, reg_carbon, reg_limiting)
            new_n_simple = mucus_simple + n_simple
            new_n_monomer = mucus_mono + n_monomer
        else
            if fix_metab
                carbon, limiting_matrix = create_metabolism(n_complex, n_simple, n_monomer, gapsize, limited_pathways=limited_pathways, seed=masterseed)
            else
                carbon, limiting_matrix = create_metabolism(n_complex, n_simple, n_monomer, gapsize, limited_pathways=limited_pathways, seed=seed)
            end
            new_n_simple = n_simple
            new_n_monomer = n_monomer
        end

        if fix_pool
            pool = create_species_pool(pool_size, n_complex, new_n_simple, new_n_monomer, carbon, seed=masterseed, limiting_matrix=limiting_matrix,
                                        gen_n_consumed=gen_n_consumed, spec_n_consumed=spec_n_consumed, consume_all_monomers=consume_all_monomers, h=h)
        else
            pool = create_species_pool(pool_size, n_complex, new_n_simple, new_n_monomer, carbon, seed=seed, limiting_matrix=limiting_matrix,
                                        gen_n_consumed=gen_n_consumed, spec_n_consumed=spec_n_consumed, consume_all_monomers=consume_all_monomers, h=h)
        end

        if fix_sample
            sample = sample_pool(pool, n_species, n_invaders, n_comms=ncomms, seed=masterseed, fix_order=fix_order)
        else
            sample = sample_pool(pool, n_species, n_invaders, n_comms=ncomms, seed=seed, fix_order=fix_order)
        end
        
        if fix_callbacks
            out, r_out = spatial_run(ncomms, sample, masterseed, phi=tradeoff_strength__number_of_splits, eta=tradeoff_strength__number_consumed, cutoff=cutoff, t_span=t_span, path=path,
                alpha_list=alpha_list, symmetrical=symmetrical, plotgraph=true, consume_all_monomers=new_n_monomer*consume_all_monomers,
                influx_amount=influx_amount, influx_time=influx_time, mucus=mucus, n_complex=n_complex, t_inv = t_inv)  
        else
            out, r_out = spatial_run(ncomms, sample, seed, phi=tradeoff_strength__number_of_splits, eta=tradeoff_strength__number_consumed, cutoff=cutoff, t_span=t_span, path=path,
                alpha_list=alpha_list, symmetrical=symmetrical, plotgraph=true, consume_all_monomers=new_n_monomer*consume_all_monomers,
                influx_amount=influx_amount, influx_time=influx_time, mucus=mucus, n_complex=n_complex, t_inv = t_inv)
        end
        ########## GETTING DATA ###########
        comm_tots, counts, present, gamma, beta, gen_tots, spec_tots = describe_full_comm_species(out)
        global t = out[1].coldata.time
        present_matrices_ts = get_present_production_matrices(present, out[1].rowdata.matrices)
        MCI = total_metabolism(present_matrices_ts)
        println("size MCI", size(MCI))
        println("vec", size(vec(MCI)))

        ########## PLOTTING ###########
        plot_species_counts(counts, t, path, false)
        plot_shannon(gamma, t, path, false)
        plot_betadiv(beta, t, path, false)
        plot_MCI(MCI, t, path, false)
        plot_consumer_types(gen_tots, spec_tots, t, path, false)
        plot_comm_totals(comm_tots, t, path, false)

        push!(total_counts, counts)
        push!(total_shannons, gamma)
        push!(total_beta, beta)
        push!(total_MCI, MCI)
        push!(total_gen, gen_tots)
        push!(total_spec, spec_tots)
        for (idx, value) in enumerate(comm_tots)
            push!(total_comm_tots[idx], value)
            println("Sum$idx:", sum(value))
        end

        gen_aggr, spec_aggr, counts_aggr, shannons_aggr, MCI_aggr = grab_ends(
            t, t_span[2], t_inv * n_invaders + 1, gen_tots, spec_tots, counts, gamma, MCI)

        for i in 1:2
            push!(gen_tr_rel[i], gen_aggr[i])
            push!(spec_tr_rel[i], spec_aggr[i])
            push!(counts_tr_rel[i], counts_aggr[i])
            push!(shannons_tr_rel[i], shannons_aggr[i])
            push!(MCI_tr_rel[i], MCI_aggr[i])
        end
    end

    plot_shannon(total_shannons, t, basepath, true, comm_id = "total")
    plot_species_counts(total_counts, t, basepath, true)
    plot_betadiv(total_beta, t, basepath, true)
    plot_MCI(total_MCI, t, basepath, true)
    plot_consumer_types(total_gen, total_spec, t, basepath, true)
    plot_comm_totals(total_comm_tots, t, basepath, true)

    names = ["gen", "spec", "counts", "shannons", "MCI"]
    metrics = [total_gen, total_spec, total_counts, total_shannons, total_MCI]

    for (name, metric) in zip(names, metrics)
        tr_rel_means_with_sd(metric, name, basepath)
    end

    @save "$basepath/gen_tr_rel.jld2" gen_tr_rel
    @save "$basepath/spec_tr_rel.jld2" spec_tr_rel
    @save "$basepath/counts_tr_rel.jld2" counts_tr_rel
    @save "$basepath/shannons_tr_rel.jld2" shannons_tr_rel
    @save "$basepath/MCI_tr_rel.jld2" MCI_tr_rel

    @save "$basepath/total_shannons.jld2" total_shannons
    @save "$basepath/total_counts.jld2" total_counts
    @save "$basepath/total_beta.jld2" total_beta
    @save "$basepath/total_MCI.jld2" total_MCI
    @save "$basepath/total_gen.jld2" total_gen
    @save "$basepath/total_spec.jld2" total_spec
    @save "$basepath/total_comm_tots.jld2" total_comm_tots
    @save "$basepath/seeds.jld2" seeds
end