function spatial_equations(u, p, t; ph_list=nothing)
    n_species = p.n_species
    n_resources = p.n_resources
    N = length(u) ÷ (n_species + n_resources)  # Number of communities

    u_mat = zeros(n_species + n_resources, N)
    u_mat[1:n_species, :] = reshape(u[1:(n_species*N)], n_species, N)
    u_mat[(n_species+1):end, :] = reshape(u[(n_species*N+1):end], n_resources, N)
    du_mat = zeros(n_species + n_resources, N)

    # First community recieves resources, the rest do not
    p_nonfirst = deepcopy(p)
    p_nonfirst.alpha = zeros(Float64, length(p.alpha))

    # FIRST COMMUNITY
    if ph_list === nothing
        du_mat[:, 1] .= equations(u_mat[:, 1], p, t)
    else
        du_mat[:, 1] .= equations(u_mat[:, 1], p, t, ph=ph_list[1])
    end
    # OTHER COMMUNITIES
    if ph_list === nothing
        for i in 2:N
            du_mat[:, i] .= equations(u_mat[:, i], p_nonfirst, t)
        end
    else
        for i in 2:N
            du_mat[:, i] .= equations(u_mat[:, i], p_nonfirst, t, ph=ph_list[i])
        end
    end 
    # -------Diffusion, Advection--------

    # MIDDLING COMMUNITIES
    for i in 2:(N-1)
        # Diffusion
        du_mat[1:n_species, i] += p.D_species .* (u_mat[1:n_species, i-1] - 2u_mat[1:n_species, i] + u_mat[1:n_species, i+1])
        du_mat[n_species+1:end, i] += p.D_resources .* (u_mat[n_species+1:end, i-1] - 2u_mat[n_species+1:end, i] + u_mat[n_species+1:end, i+1])

        # Advection
        du_mat[1:n_species, i] += p.A_species .* u_mat[1:n_species, i-1] - p.A_species .* u_mat[1:n_species, i]
        du_mat[n_species+1:end, i] += p.A_resources .* u_mat[n_species+1:end, i-1] - p.A_resources .* u_mat[n_species+1:end, i]
    end

    # FIRST COMMUNITY
    # Diffusion
    du_mat[1:n_species, 1] -= p.A_species .* u_mat[1:n_species, 1]
    du_mat[n_species+1:end, 1] -= p.A_resources .* u_mat[n_species+1:end, 1]

    # Advection
    du_mat[1:n_species, 1] += p.D_species .* (u_mat[1:n_species, 2] - u_mat[1:n_species, 1])
    du_mat[n_species+1:end, 1] += p.D_resources .* (u_mat[n_species+1:end, 2] - u_mat[n_species+1:end, 1])

    # LAST COMMUNITY
    # Diffusion
    du_mat[1:n_species, N] += p.D_species .* (u_mat[1:n_species, N-1] - u_mat[1:n_species, N])
    du_mat[n_species+1:end, N] += p.D_resources .* (u_mat[n_species+1:end, N-1] - u_mat[n_species+1:end, N])

    # Advection
    du_mat[1:n_species, N] += p.A_species .* u_mat[1:n_species, N-1] - p.A_species .* u_mat[1:n_species, N]
    du_mat[n_species+1:end, N] += p.A_resources .* u_mat[n_species+1:end, N-1] - p.A_resources .* u_mat[n_species+1:end, N]

    # ----------------------------------

    # Flatten du_mat into du
    du = zeros((n_species+n_resources)*N)
    du[1:(n_species*N)] .= reshape(du_mat[1:n_species, :], n_species*N)
    du[(n_species*N+1):end] .= reshape(du_mat[(n_species+1):end, :], n_resources*N)
    return du
end