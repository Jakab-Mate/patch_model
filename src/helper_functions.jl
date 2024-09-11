function non_zero_pos(A)
    indices = findall(x -> x != 0, A)
    tuple_indices = Tuple.(indices)
    return tuple_indices
end

function normalize(A; h::Real=0)
    if isempty(A)
        throw(DomainError("Cannot normalize an empty array."))
    else
        if sum(A) == 0
            throw(DomainError("Cannot normalize an array that sums to zero."))
        end
    end

    if h == 0
        return A ./ sum(A)
    else
        return A ./ sum(A .^ h)
    end
           
end

function sample_reaction_indices(rng, D, number_of_reactions)
    total_different_reactions = length(non_zero_pos(D))
    if total_different_reactions < number_of_reactions
        throw(DomainError("Number of sampled reactions ($number_of_reactions) should be smaller than the total number of reactions ($total_different_reactions)"))
    end

    helper_D = copy(D)
    chosen_indices = Array{Tuple}(undef, number_of_reactions)
    for i in 1:number_of_reactions
        choice = rand(rng, non_zero_pos(helper_D))
        chosen_indices[i] = choice
        #excluding production and consumption of the same resource
        #helper_D[:, choice[1]] .= 0
        #helper_D[choice[2], :] .= 0
        helper_D[choice[1], choice[2]] = 0
    end

    return chosen_indices
end

function checks_before_run(D, W_ba, tau, alpha)
    if isnothing(D)
        if !isnothing(W_ba)
            @warn "WARNING: Supplied energy yield matrix (W_ba) but no stoichiometric matrix (D). Creating D matrix of same size"
            D, W_ba = create_metabolism(n_resources=length(W_ba[1, :]))
        else
            D, W_ba = create_metabolism()
        end
    else
        if isnothing(W_ba)
            D, W_ba = create_metabolism(n_resources=length(D[1, :]))
            @warn "WARNING: Supplied stoichiometric matrix (D) but no energy yield matrix (W_ba). Creating W_ba matrix of same size"
        end
    end

    n_resources = size(D, 1)

    D_row, D_col = size(D)
    if (n_resources != D_row) || (n_resources != D_col)
        throw(DomainError("Number of resources does not match the size of the stoichiometric matrix (D) \n 
        number of resources is $n_resources, sizes of D are ($D_row, $D_col)"))
    end

    W_row, W_col = size(W_ba)
    if (n_resources != W_row) || (n_resources != W_col)
        throw(DomainError("Number of resources does not match the size of the energy yield matrix (W_ba) \n 
        number of resources is $n_resources, sizes of D are ($W_row, $W_col)"))
    end

    if isnothing(tau)
        tau = ones(Float64, n_resources)
    end

    if isnothing(alpha)
        alpha = vcat([100.0], zeros(Float64, n_resources-1))
    end

    if length(tau) != n_resources || length(alpha) != n_resources
        throw(DomainError("Length of dilution terms (tau) or resource availabilities (alpha) does not match the number of resources \n This can happen when tau or alpha is supplied, but the stoichiometric matrix and energy yield matrix are not. \n"))
    end

    return n_resources, D, W_ba, tau, alpha
end