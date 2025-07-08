function aggregate_by_rowdata(se, row, assay_name)
    row_data = se.rowdata[!, Symbol(row)]
    row_set = Set(row_data)
    ts_data = Number.(assay(se, assay_name))'
    aggregated_data = Dict()
    for unique_val in row_set
        indices = findall(x -> x == unique_val, row_data)
        summed_vector = sum(ts_data[:, indices], dims=2)
        aggregated_data[unique_val] = summed_vector
    end

    return aggregated_data
end

function my_mean(arr)
    return sum(arr) / length(arr)
end

function my_std(arr)
    m = my_mean(arr)
    return sqrt(sum((x - m)^2 for x in arr) / (length(arr) - 1))
end

function plot_se(se::SummarizedExperiment, assay_name::String, path::String; comm=false)
    labels = reshape(se.rowdata.name, (1, length(se.rowdata.name)))
    time = se.coldata.time
    data = Number.(assay(se, assay_name))

    p = plot(time, Number.(assay(se, assay_name))',
        label = labels, legend_position = :outerleft, xaxis = "Time", yaxis = "Abundance")
    
    if !isinteger(comm)
        savefig(p, joinpath(path, "strain_abundances.png"))
    else
        savefig(p, joinpath(path, "comm$comm" * "strain_abundances.png"))
    end

    for attr in ["n_reactions"]
        #if haskey(se.rowdata, Symbol(attr))
            aggreg_data = aggregate_by_rowdata(se, attr, assay_name)
            aggreg_array = collect(keys(aggreg_data))
            labels = reshape(aggreg_array, (1, length(aggreg_array)))
            p = plot(se.coldata.time, [x for x in values(aggreg_data)],
                label = labels, legend_position = :outerleft, xaxis = "Time", yaxis = "Abundance")
            if !isinteger(comm)
                savefig(p, joinpath(path, "$attr _abundances.png"))
            else
                savefig(p, joinpath(path, "comm$comm" * "$attr _abundances.png"))
            end
        #end
    end
end

function plot_community_graph(consums::Vector{Array{Float64}}, prods, consumption_fluxes, production_fluxes, species_abund,
     path::String; consume_all_monomers::Int64=0)

    n_communities = length(consumption_fluxes)
    println("Number of communities: ", n_communities)
    n_species = length(prods)
    println("Number of species: ", n_species)
    if n_species == 0
        println("No species survived")
        return
    end
    n_resources = size(prods[1], 1)

    if consume_all_monomers > 0
        n_resources -= consume_all_monomers
        # for comm in 1:n_communities
        #     consumption_fluxes[comm] = [x[1:end-consume_all_monomers] for x in consumption_fluxes[comm]]
        #     production_fluxes[comm] = [x[1:end-consume_all_monomers] for x in production_fluxes[comm]]
        # end
        ### consums and prods should be updated as well, but it might all be unnecessary, because decereasing n_resources ensures we never get to the monomer indices
    end

    for comm in 1:n_communities
        g = SimpleDiGraph()

        for i in 1:(n_species + n_resources)
            add_vertex!(g)
        end

        edge_weights = Dict{Tuple{Int, Int}, Float64}()

        for species_idx in 1:n_species
            for col in 1:n_resources
                if consums[species_idx][col] != 0.0
                    add_edge!(g, n_species + col, species_idx)  # Red arrow
                    edge_weights[(n_species + col, species_idx)] = consumption_fluxes[comm][species_idx][col]
                end
                for row in 1:n_resources
                    if prods[species_idx][row, col] != 0.0
                        add_edge!(g, species_idx, n_species + row)  # Green arrow
                        edge_weights[(species_idx, n_species + row)] = production_fluxes[comm][species_idx][row]
                    end
                end
            end
        end

        species_labels = ["S$i" for i in 1:n_species]
        resource_labels = ["R$i" for i in 1:n_resources]
        labels = vcat(species_labels, resource_labels)
        node_colors = vcat(fill("lightblue", n_species), fill("lightgreen", n_resources))

        edgecolors = []
        edgestrokes = []
        for edge in edges(g)
            src, dst = edge.src, edge.dst
            if src <= n_species && dst > n_species
                push!(edgecolors, "green")
            else
                push!(edgecolors, "red")
            end
            # Use edge weights for stroke thickness
            weight = edge_weights[(src, dst)]
            if weight > 0
                weight = log(weight)
            else
                weight = 0.0
            end

            push!(edgestrokes, 0.5 + 1.0 * weight)
        end

        abundances = species_abund[comm]
        log_abundances = []
        for abundance in abundances
            if abundance > 0
                push!(log_abundances, log(abundance))
            else
                push!(log_abundances, 0.0)
            end
        end
        resource_weight = maximum(log_abundances)
        node_sizes = vcat(log_abundances, fill(resource_weight, n_resources))

        plot_g = gplot(
            g,
            layout=circular_layout,
            nodelabel=labels,
            nodefillc=node_colors,
            nodesize=node_sizes,
            edgelinewidth=edgestrokes,
            edgestrokec=edgecolors
        )

        background = compose(plot_g, 
                        rectangle(-1.2, -1.2, 2.4, 2.4), 
                        fill("white"))
        
        final_plot = compose(background, plot_g)
        println("Saving community graph $comm at location $path")
        draw(PNG(joinpath(path, "community_graph$comm.png"), 16cm, 16cm), final_plot)
    end
end


function plot_community(end_result)
    # create graph using LightGraphs.jl
end

function plot_shannon(shannons, t, path, multiple; comm_id = "")
    println("Size shannons: ", size(shannons))
    if multiple
        data = reduce(hcat, shannons)
        means = [sum(row) / length(row) for row in eachrow(data)]
        stds = [sqrt(sum((x - mean)^2 for x in row) / length(row)) for (row, mean) in zip(eachrow(data), means)]

        fig = plot(t, means, ribbon=stds, fillalpha=0.2, xlabel="Time", ylabel="Gamma diversity", linewidth=3)
    else
        fig = plot(t, shannons, xlabel="Time", ylabel="Gamma diversity", linewidth=3)
    end
    savefig(fig, "$path/shannon_diversity_$comm_id.png")
end

function plot_species_counts(counts, t, path, multiple)
    println("Size counts: ", size(counts))
    if multiple
        data = reduce(hcat, counts)
        means = [sum(row) / length(row) for row in eachrow(data)]
        stds = [sqrt(sum((x - mean)^2 for x in row) / length(row)) for (row, mean) in zip(eachrow(data), means)]
        fig = plot(t, means, ribbon=stds, fillalpha=0.2, xlabel="Time", ylabel="Number of different species", linewidth=3)
    
    else
        fig = plot(t, counts, linewidth=3, xlabel="Time", ylabel="Number of different species")
    
    end
    savefig(fig, "$path/species_counts")
end

function plot_betadiv(betadivs, t, path, multiple)
    if multiple
        data = reduce(hcat, betadivs)
        means = [sum(row) / length(row) for row in eachrow(data)]
        stds = [sqrt(sum((x - mean)^2 for x in row) / length(row)) for (row, mean) in zip(eachrow(data), means)]
        fig = plot(t, means, ribbon=stds, fillalpha=0.2, xlabel="Time", ylabel="Beta diversity", linewidth=3, ylim=(-1, 3))
    else
        fig = plot(t, betadivs, linewidth=3, xlabel="Time", ylabel="Beta diversity", ylim=(-1, 3))
    end
    savefig(fig, "$path/beta_diversity")
end

function plot_MCI(MCI_ts, t, path, multiple)
    if multiple
        data = reduce(hcat, MCI_ts)
        means = [sum(row) / length(row) for row in eachrow(data)]
        stds = [sqrt(sum((x - mean)^2 for x in row) / length(row)) for (row, mean) in zip(eachrow(data), means)]
        fig = plot(t, means, ribbon=stds, fillalpha=0.2, xlabel="Time", ylabel="MCI", linewidth=3)
    else
        fig = plot(t, MCI_ts, linewidth=3, xlabel="Time", ylabel="Number of diff. reactions")
    end
    savefig(fig, "$path/MCI")
end

function plot_consumer_types(gen_ts, spec_ts, t, path, multiple)
    if multiple
        gen_data = reduce(hcat, gen_ts)
        gen_means = [sum(row) / length(row) for row in eachrow(gen_data)]
        #gen_stds = [sqrt(sum((x - mean)^2 for x in row) / length(row)) for (row, mean) in zip(eachrow(gen_data), gen_means)]
        spec_data = reduce(hcat, spec_ts)
        spec_means = [sum(row) / length(row) for row in eachrow(spec_data)]
        #spec_stds = [sqrt(sum((x - mean)^2 for x in row) / length(row)) for (row, mean) in zip(eachrow(spec_data), spec_means)]
        fig = plot(t, [gen_means, spec_means], xlabel="Time", ylabel="Total abundance",
            label=["Primary consumers" "Secondary consumers"], linewidth=3)
    else
        fig = plot(t, [gen_ts, spec_ts], linewidth=3, xlabel="Time", ylabel="Total abundance",
            label=["Primary consumers" "Secondary consumers"])
    end
    savefig(fig, "$path/consumer_types")
end

function plot_comm_totals(comm_tot_abundances_ts, t, path, multiple)
    labels_0 = ["comm $i" for i in 1:length(comm_tot_abundances_ts)]
    labels = reshape(labels_0, (1, length(labels_0)))
    if multiple
        means = []
        stds = []
        for comm in comm_tot_abundances_ts
            println("Size comm: ", size(comm))
            data = reduce(hcat, comm)
            println("Size data: ", size(data))
            curr_means = [sum(row) / length(row) for row in eachrow(data)]
            push!(means, curr_means)
            #curr_stds = [sqrt(sum((x - mean)^2 for x in row) / length(row)) for (row, mean) in zip(eachrow(data), curr_means)]
            #push!(stds, curr_stds)
        end
        println("Size means: ", size(means))
        println("Size means[1]: ", size(means[1]))
        println("Size means[2]: ", size(means[2]))
        println("Size means[3]: ", size(means[3]))
        fig = plot(t, means, xlabel="Time", ylabel="Total abundance", linewidth=3, label=labels)
    else
        fig = plot(t, comm_tot_abundances_ts, linewidth=3, xlabel="Time", ylabel="Total abundance", label=labels)
    end
    savefig(fig, "$path/comm_totals")
end

function tr_rel_means_with_sd(data, name, path)
    tr = data[1]
    rel = data[2]

    means = [my_mean(tr), my_mean(rel)]
    stds = [my_std(tr), my_std(rel)]

    x = [1, 2]

    plot(x, means, yerror=stds, seriestype=:scatter, markershape=:circle, color=:blue)
    plot!(x, means, linewidth=2, linestyle=:dash, label="")

    xticks!(x, ["Transient", "Relaxation"])  # Set x-axis labels
    ylabel!(name)
    title!(uppercasefirst(name) * " averages at transient and relaxation ends")

    savefig(joinpath(path, "tr_rel_$name.png"))
end