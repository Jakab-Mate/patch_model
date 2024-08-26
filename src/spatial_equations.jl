function spatial_equations(u, p, t)
    n_species = p.n_species
    n_resources = p.n_resources
    N = length(u) ÷ (n_species + n_resources)  # Number of communities

    u_mat = zeros(n_species + n_resources, N)
    u_mat[1:n_species, :] = reshape(u[1:(n_species*N)], n_species, N)
    u_mat[(n_species+1):end, :] = reshape(u[(n_species*N+1):end], n_resources, N)
    du_mat = zeros(n_species + n_resources, N)

    for i in 1:N
        du_mat[:, i] .= equations(u_mat[:, i], p, t)
    end

    for i in 2:(N-1)
        du_mat[1:n_species, i] += p.D_species .* (u_mat[1:n_species, i-1] - 2u_mat[1:n_species, i] + u_mat[1:n_species, i+1])
        du_mat[n_species+1:end, i] += p.D_resources .* (u_mat[n_species+1:end, i-1] - 2u_mat[n_species+1:end, i] + u_mat[n_species+1:end, i+1])
    end

    # Add unidirectional flow

    du_mat[1:n_species, 1] += p.D_species .* (u_mat[1:n_species, 2] - u_mat[1:n_species, 1])
    du_mat[n_species+1:end, 1] += p.D_resources .* (u_mat[n_species+1:end, 2] - u_mat[n_species+1:end, 1])

    du_mat[1:n_species, N] += p.D_species .* (u_mat[1:n_species, N-1] - u_mat[1:n_species, N])
    du_mat[n_species+1:end, N] += p.D_resources .* (u_mat[n_species+1:end, N-1] - u_mat[n_species+1:end, N])

    du = zeros((n_species+n_resources)*N)
    du[1:(n_species*N)] .= reshape(du_mat[1:n_species, :], n_species*N)
    du[(n_species*N+1):end] .= reshape(du_mat[(n_species+1):end, :], n_resources*N)
    return du
end