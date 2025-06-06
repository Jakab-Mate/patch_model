function repeat_params(n_repeats, basepath; 
    n_complex=3, #ok
    n_simple=10, #ok
    n_monomer=5, #ok
    gapsize=15, #ok
    n_resources = n_complex + n_simple + n_monomer, #?
    n_species = 1, # set to 0? ok
    n_invaders = 40,  # 220 ok
    pool_size = 500, #ok
    ncomms = 3, #ok
    t_span = (0, 60), #ok
    cutoff = 0.0001, #ok
    tradeoff_strength__number_consumed = 0.1, #ok
    tradeoff_strength__number_of_splits = 0.1, #ok
    gen_n_consumed = 4, #ok
    spec_n_consumed = 1, #ok
    symmetrical = true, #ok
    limited_pathways = true, #ok
    consume_all_monomers = true, #ok
    h = 1.0, #ok
    influx_amount = 0.0, #ok
    influx_time = 1.0, #ok
    mucus = true, #ok
    mucus_disjunct = true, #ok
    mucus_simple = 5, #ok
    mucus_mono = 4, #ok
    mucus_gap = 10, #ok
    t_inv = 1.0, #ok
    alpha_high = 100.0, #ok
    set_seeds = nothing, #ok
    fix_metab = false, #ok
    fix_pool = false, #ok
    fix_sample = false, #ok
    fix_callbacks = false,
    shuffle = false) #ok

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
            sample = sample_pool(pool, n_species, n_invaders, n_comms=ncomms, seed=masterseed, shuffle=shuffle)
        else
            sample = sample_pool(pool, n_species, n_invaders, n_comms=ncomms, seed=seed, shuffle=shuffle)
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