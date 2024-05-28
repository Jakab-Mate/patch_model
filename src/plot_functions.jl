using MicrobiomeAnalysis, Plots

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

function plot_MCI(se::SummarizedExperiment, assay_name::String, path::String)
    ts_data = Number.(assay(se, assay_name))'
    n_react_array = se.rowdata[!, Symbol("n_reactions")]
    present_accumulator = zeros(Int64, length(se.coldata.time))
    react_accumulator = zeros(Int64, length(se.coldata.time))
    for i in eachindex(n_react_array)
        arr = ts_data[:, i]
        index = findfirst(x -> !isdigit(x), n_react_array[i])
        num = parse(Int, n_react_array[i][1:index-1])
        arr[arr .> 0.0001] .= num
        arr[arr .<= 0.0001] .= 0
        react_accumulator .+= arr
        arr[arr .> 0.0001] .= 1
        present_accumulator .+= arr
    end
    
    arr_names = ["Cumulative", "Per Capita"]
    plotter = hcat(react_accumulator, react_accumulator ./ present_accumulator)
    labels = reshape(arr_names, (1, 2))

    p = Plots.plot(se.coldata.time, plotter,
        label = labels, legend_position = :outerleft, xaxis = "Time", yaxis = "MCI")
    savefig(p, joinpath(path, "MCI.png"))
    println("plot was saved at $path")
end

function plot_se(se::SummarizedExperiment, assay_name::String, path::String)
    labels = reshape(se.rowdata.name, (1, length(se.rowdata.name)))

    p = Plots.plot(se.coldata.time, Number.(assay(se, assay_name))',
        label = labels, legend_position = :outerleft, xaxis = "Time", yaxis = "Abundance")

    savefig(p, joinpath(path, "strain_abundances.png"))

    for attr in ["family", "n_reactions", "n_splits"]
        aggreg_data = aggregate_by_rowdata(se, attr, assay_name)
        aggreg_array = collect(keys(aggreg_data))
        labels = reshape(aggreg_array, (1, length(aggreg_array)))
        p = Plots.plot(se.coldata.time, [x for x in values(aggreg_data)],
            label = labels, legend_position = :outerleft, xaxis = "Time", yaxis = "Abundance")
        savefig(p, joinpath(path, "$attr _abundances.png"))
    end
end


function plot_community(end_result)
    # create graph using LightGraphs.jl
end