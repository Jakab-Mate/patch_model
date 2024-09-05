function affect!(integrator, t_inv, start_time, n_invaders, n_species, cutoff)
    integrator.p.present_species = findall(x -> x > cutoff, integrator.u[1:(n_species+n_invaders)])
    if integrator.t >= start_time
        idx = Int(n_species + floor((integrator.t - start_time) / t_inv)) + 1
        if idx <= n_invaders + n_species
            println("invader ", idx-n_species, " added at time ", integrator.t)
            integrator.u[idx] += 10
        end
    end
end

function create_callbacks(t_inv, start_time, n_invaders, n_species, cutoff)
    cb = PeriodicCallback(
        integrator -> affect!(integrator, t_inv, start_time, n_invaders, n_species, cutoff), t_inv
    )
    return cb
end

# Should cutoff be divided by number of communities?

function spatial_affect!(integrator, t_inv, start_time, n_invaders, n_species, cutoff, n_comms)
    if integrator.t >= start_time
        idx = Int(n_species + floor((integrator.t - start_time) / t_inv)) + 1
        if idx <= n_invaders + n_species
            println("invader ", idx-n_species, " added at time ", integrator.t)
            integrator.u[idx] += 10 # Adding invader only to first community
        end
    end

    present_species_set = Set{Int}()
    for i in 1:n_comms
        present_in_current = findall(x -> x > cutoff, integrator.u[(i-1)*(n_species+n_invaders)+1:i*(n_species+n_invaders)])
        present_species_set = union(present_species_set, present_in_current)
    end
    integrator.p.present_species = collect(present_species_set)
end

function create_spatial_callbacks(t_inv, start_time, n_invaders, n_species, cutoff, n_comms)
    cb = PeriodicCallback(
        integrator -> spatial_affect!(integrator, t_inv, start_time, n_invaders, n_species, cutoff, n_comms), t_inv
    )
    return cb
end

