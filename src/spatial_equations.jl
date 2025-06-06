@views function spatial_equations!(du, u, p, t; N=nothing, u_mat=nothing, du_mat=nothing, ph_list=nothing, alpha_list=nothing, symmetrical=false)
    u_mat[1:p.n_species, :] = reshape(u[1:(p.n_species*N)], p.n_species, N)
    u_mat[(p.n_species+1):end, :] = reshape(u[(p.n_species*N+1):end], p.n_resources, N)
    
    if ph_list === nothing
        for i in 1:N
            du_mat[:, i] = threaded_equations(u_mat[:, i], p, t, alpha=alpha_list[i])
        end
    else
        for i in 1:N
            du_mat[:, i] = threaded_equations(u_mat[:, i], p, t, ph=ph_list[i], alpha=alpha_list[i])
        end
    end 

    # -------Diffusion, Advection--------
    if symmetrical
        # Only diffusion
        # MIDDLING COMMUNITIES
        for i in 2:N-1
            du_mat[1:p.n_species, i] += p.D_species .* (u_mat[1:p.n_species, i-1] - 2 * u_mat[1:p.n_species, i] + u_mat[1:p.n_species, i+1])
            du_mat[p.n_species+1:end, i] += p.D_resources .* (u_mat[p.n_species+1:end, i-1] - 2 * u_mat[p.n_species+1:end, i] + u_mat[p.n_species+1:end, i+1])
        end

        # FIRST COMMUNITY
        du_mat[1:p.n_species, 1] += p.D_species .* (u_mat[1:p.n_species, 2] - 2 * u_mat[1:p.n_species, 1] + u_mat[1:p.n_species, N])
        du_mat[p.n_species+1:end, 1] += p.D_resources .* (u_mat[p.n_species+1:end, 2] - 2 * u_mat[p.n_species+1:end, 1] + u_mat[p.n_species+1:end, N])

        # LAST COMMUNITY
        du_mat[1:p.n_species, N] += p.D_species .* (u_mat[1:p.n_species, 1] - 2 * u_mat[1:p.n_species, N] + u_mat[1:p.n_species, N-1])
        du_mat[p.n_species+1:end, N] += p.D_resources .* (u_mat[p.n_species+1:end, 1] - 2 * u_mat[p.n_species+1:end, N] + u_mat[p.n_species+1:end, N-1])
    else
        # MIDDLING COMMUNITIES
        for i in 2:(N-1)
            # Diffusion
            du_mat[1:p.n_species, i] += p.D_species .* (u_mat[1:p.n_species, i-1] - 2 * u_mat[1:p.n_species, i] + u_mat[1:p.n_species, i+1])
            du_mat[p.n_species+1:end, i] += p.D_resources .* (u_mat[p.n_species+1:end, i-1] - 2 * u_mat[p.n_species+1:end, i] + u_mat[p.n_species+1:end, i+1])

            # Advection
            du_mat[1:p.n_species, i] += p.A_species .* u_mat[1:p.n_species, i-1] - p.A_species .* u_mat[1:p.n_species, i]
            du_mat[p.n_species+1:end, i] += p.A_resources .* u_mat[p.n_species+1:end, i-1] - p.A_resources .* u_mat[p.n_species+1:end, i]
        end

        # FIRST COMMUNITY
        # Advection
        du_mat[1:p.n_species, 1] -= p.A_species .* u_mat[1:p.n_species, 1]
        du_mat[p.n_species+1:end, 1] -= p.A_resources .* u_mat[p.n_species+1:end, 1]

        # Diffusion
        du_mat[1:p.n_species, 1] += p.D_species .* (u_mat[1:p.n_species, 2] - u_mat[1:p.n_species, 1])
        du_mat[p.n_species+1:end, 1] += p.D_resources .* (u_mat[p.n_species+1:end, 2] - u_mat[p.n_species+1:end, 1])

        # LAST COMMUNITY
        # Diffusion
        du_mat[1:p.n_species, N] += p.D_species .* (u_mat[1:p.n_species, N-1] - u_mat[1:p.n_species, N])
        du_mat[p.n_species+1:end, N] += p.D_resources .* (u_mat[p.n_species+1:end, N-1] - u_mat[p.n_species+1:end, N])

        # Advection
        du_mat[1:p.n_species, N] += p.A_species .* u_mat[1:p.n_species, N-1] - p.A_species .* u_mat[1:p.n_species, N]
        du_mat[p.n_species+1:end, N] += p.A_resources .* u_mat[p.n_species+1:end, N-1] - p.A_resources .* u_mat[p.n_species+1:end, N]
    end
    # ----------------------------------

    # Flatten du_mat into du
    # du = zeros((p.n_species+p.n_resources)*N)
    du[1:(p.n_species*N)] .= reshape(du_mat[1:p.n_species, :], p.n_species*N)
    du[(p.n_species*N+1):end] .= reshape(du_mat[(p.n_species+1):end, :], p.n_resources*N)
    nothing
end

@views function essential_spatial_equations!(du, u, p, t; N=nothing, u_mat=nothing, du_mat=nothing, ph_list=nothing, alpha_list=nothing, symmetrical=false)
    u_mat[1:p.n_species, :] = reshape(u[1:(p.n_species*N)], p.n_species, N)
    u_mat[(p.n_species+1):end, :] = reshape(u[(p.n_species*N+1):end], p.n_resources+1, N)
    
    if ph_list === nothing
        for i in 1:N
            du_mat[:, i] = essential_R_equations(u_mat[:, i], p, t, alpha=alpha_list[i])
        end
    else
        for i in 1:N
            du_mat[:, i] = essential_R_equations(u_mat[:, i], p, t, ph=ph_list[i], alpha=alpha_list[i])
        end
    end 

    # -------Diffusion, Advection--------

    # MIDDLING COMMUNITIES
    for i in 2:(N-1)
        # Diffusion
        du_mat[1:p.n_species, i] += p.D_species .* (u_mat[1:p.n_species, i-1] - 2 * u_mat[1:p.n_species, i] + u_mat[1:p.n_species, i+1])
        du_mat[p.n_species+1:end, i] += p.D_resources .* (u_mat[p.n_species+1:end, i-1] - 2 * u_mat[p.n_species+1:end, i] + u_mat[p.n_species+1:end, i+1])

        # Advection
        du_mat[1:p.n_species, i] += p.A_species .* u_mat[1:p.n_species, i-1] - p.A_species .* u_mat[1:p.n_species, i]
        du_mat[p.n_species+1:end, i] += p.A_resources .* u_mat[p.n_species+1:end, i-1] - p.A_resources .* u_mat[p.n_species+1:end, i]
    end

    # FIRST COMMUNITY
    # Diffusion
    du_mat[1:p.n_species, 1] -= p.A_species .* u_mat[1:p.n_species, 1]
    du_mat[p.n_species+1:end, 1] -= p.A_resources .* u_mat[p.n_species+1:end, 1]

    # Advection
    du_mat[1:p.n_species, 1] += p.D_species .* (u_mat[1:p.n_species, 2] - u_mat[1:p.n_species, 1])
    du_mat[p.n_species+1:end, 1] += p.D_resources .* (u_mat[p.n_species+1:end, 2] - u_mat[p.n_species+1:end, 1])

    # LAST COMMUNITY
    # Diffusion
    du_mat[1:p.n_species, N] += p.D_species .* (u_mat[1:p.n_species, N-1] - u_mat[1:p.n_species, N])
    du_mat[p.n_species+1:end, N] += p.D_resources .* (u_mat[p.n_species+1:end, N-1] - u_mat[p.n_species+1:end, N])

    # Advection
    du_mat[1:p.n_species, N] += p.A_species .* u_mat[1:p.n_species, N-1] - p.A_species .* u_mat[1:p.n_species, N]
    du_mat[p.n_species+1:end, N] += p.A_resources .* u_mat[p.n_species+1:end, N-1] - p.A_resources .* u_mat[p.n_species+1:end, N]

    # ----------------------------------

    # Flatten du_mat into du
    # du = zeros((p.n_species+p.n_resources)*N)
    du[1:(p.n_species*N)] .= reshape(du_mat[1:p.n_species, :], p.n_species*N)
    du[(p.n_species*N+1):end] .= reshape(du_mat[(p.n_species+1):end, :], (p.n_resources+1)*N)
    nothing
end


@views function saving_spatial_equations!(du, u, p, t; N=nothing, u_mat=nothing, du_mat=nothing, ph_list=nothing, alpha_list=nothing, symmetrical=false)
    u_mat[1:p.n_species, :] = reshape(u[1:(p.n_species*N)], p.n_species, N)
    u_mat[(p.n_species+1):end, :] = reshape(u[(p.n_species*N+1):end], p.n_resources, N)
    
    if ph_list === nothing
        for i in 1:N
            p.consumption_fluxes[1] = []
            p.production_fluxes[1] = []
            du_mat[:, i] = saving_equations(u_mat[:, i], p, t, alpha=alpha_list[i])
            p.consumption_fluxes[i+1] = p.consumption_fluxes[1]
            p.production_fluxes[i+1] = p.production_fluxes[1]
        end
    else
        for i in 1:N
            p.consumption_fluxes[1] = []
            p.production_fluxes[1] = []
            du_mat[:, i] = saving_equations(u_mat[:, i], p, t, ph=ph_list[i], alpha=alpha_list[i])
            p.consumption_fluxes[N+1] = p.consumption_fluxes[1]
            p.production_fluxes[N+1] = p.production_fluxes[1]
        end
    end 

    # -------Diffusion, Advection--------
    if symmetrical
        # Only diffusion
        # MIDDLING COMMUNITIES
        for i in 2:N-1
            du_mat[1:p.n_species, i] += p.D_species .* (u_mat[1:p.n_species, i-1] - 2 * u_mat[1:p.n_species, i] + u_mat[1:p.n_species, i+1])
            du_mat[p.n_species+1:end, i] += p.D_resources .* (u_mat[p.n_species+1:end, i-1] - 2 * u_mat[p.n_species+1:end, i] + u_mat[p.n_species+1:end, i+1])
        end

        # FIRST COMMUNITY
        du_mat[1:p.n_species, 1] += p.D_species .* (u_mat[1:p.n_species, 2] - 2 * u_mat[1:p.n_species, 1] + u_mat[1:p.n_species, N])
        du_mat[p.n_species+1:end, 1] += p.D_resources .* (u_mat[p.n_species+1:end, 2] - 2 * u_mat[p.n_species+1:end, 1] + u_mat[p.n_species+1:end, N])

        # LAST COMMUNITY
        du_mat[1:p.n_species, N] += p.D_species .* (u_mat[1:p.n_species, 1] - 2 * u_mat[1:p.n_species, N] + u_mat[1:p.n_species, N-1])
        du_mat[p.n_species+1:end, N] += p.D_resources .* (u_mat[p.n_species+1:end, 1] - 2 * u_mat[p.n_species+1:end, N] + u_mat[p.n_species+1:end, N-1])
    else
        # MIDDLING COMMUNITIES
        for i in 2:(N-1)
            # Diffusion
            du_mat[1:p.n_species, i] += p.D_species .* (u_mat[1:p.n_species, i-1] - 2 * u_mat[1:p.n_species, i] + u_mat[1:p.n_species, i+1])
            du_mat[p.n_species+1:end, i] += p.D_resources .* (u_mat[p.n_species+1:end, i-1] - 2 * u_mat[p.n_species+1:end, i] + u_mat[p.n_species+1:end, i+1])

            # Advection
            du_mat[1:p.n_species, i] += p.A_species .* u_mat[1:p.n_species, i-1] - p.A_species .* u_mat[1:p.n_species, i]
            du_mat[p.n_species+1:end, i] += p.A_resources .* u_mat[p.n_species+1:end, i-1] - p.A_resources .* u_mat[p.n_species+1:end, i]
        end

        # FIRST COMMUNITY
        # Advection
        du_mat[1:p.n_species, 1] -= p.A_species .* u_mat[1:p.n_species, 1]
        du_mat[p.n_species+1:end, 1] -= p.A_resources .* u_mat[p.n_species+1:end, 1]

        # Diffusion
        du_mat[1:p.n_species, 1] += p.D_species .* (u_mat[1:p.n_species, 2] - u_mat[1:p.n_species, 1])
        du_mat[p.n_species+1:end, 1] += p.D_resources .* (u_mat[p.n_species+1:end, 2] - u_mat[p.n_species+1:end, 1])

        # LAST COMMUNITY
        # Diffusion
        du_mat[1:p.n_species, N] += p.D_species .* (u_mat[1:p.n_species, N-1] - u_mat[1:p.n_species, N])
        du_mat[p.n_species+1:end, N] += p.D_resources .* (u_mat[p.n_species+1:end, N-1] - u_mat[p.n_species+1:end, N])

        # Advection
        du_mat[1:p.n_species, N] += p.A_species .* u_mat[1:p.n_species, N-1] - p.A_species .* u_mat[1:p.n_species, N]
        du_mat[p.n_species+1:end, N] += p.A_resources .* u_mat[p.n_species+1:end, N-1] - p.A_resources .* u_mat[p.n_species+1:end, N]
    end
    # ----------------------------------

    # Flatten du_mat into du
    # du = zeros((p.n_species+p.n_resources)*N)
    du[1:(p.n_species*N)] .= reshape(du_mat[1:p.n_species, :], p.n_species*N)
    du[(p.n_species*N+1):end] .= reshape(du_mat[(p.n_species+1):end, :], p.n_resources*N)
    nothing
end