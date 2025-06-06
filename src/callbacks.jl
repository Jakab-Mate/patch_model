"""
Create a callback function that applies `affect!` at regular intervals.
"""
function create_callbacks(t_inv, start_time, n_invaders, n_species, cutoff)
    cb = PeriodicCallback(
        integrator -> affect!(integrator, t_inv, start_time, n_invaders, n_species, cutoff), t_inv
    )
    return cb
end


"""
Remove any species with abundance below `cutoff` from the list of present species.
If `start_time` has passed and the time since the last invader addition is greater than `t_inv`,
add an invader to a random community, by setting its previously 0 abundance to 10.
"""
function spatial_affect!(integrator, t_inv, start_time, n_invaders, n_species, cutoff, n_comms, rng)
    if integrator.t >= start_time
        comm_to_add_to = rand(rng, 1:n_comms)
        idx = Int(n_species + floor((integrator.t - start_time) / t_inv)) + 1 + (comm_to_add_to-1)*(n_species+n_invaders)
        if idx <= comm_to_add_to * (n_invaders + n_species)
            println("invader ", idx - n_species - (comm_to_add_to-1) * (n_species + n_invaders), " added at time ", integrator.t)
            integrator.u[idx] += 10 
        end
    end

    present_species_set = Set{Int}()
    for i in 1:n_comms
        present_in_current = findall(x -> x > cutoff, integrator.u[(i-1)*(n_species+n_invaders)+1:i*(n_species+n_invaders)])
        present_species_set = union(present_species_set, present_in_current)
    end
    integrator.p.present_species = collect(present_species_set)
    println("present species: ", integrator.p.present_species)
end

"""
Add `influx_amount` of a random complex resource to each community. 
"""
function resource_influx!(integrator, n_species, n_invaders, n_comms, influx_amount, n_complex, rng)
    # TODO: Exclude mucus from the influx
    # TODO: Should the complex resource be random?
    for comm in 1:n_comms
        base_idx = comm * (n_species + n_invaders)
        complex_idx = rand(rng, 1:n_complex)
        idx = base_idx + complex_idx
        integrator.u[idx] += influx_amount
    end
end

"""
Create a callback function that applies `spatial_affect!` at `t_inv` intervals.
Create a second callback function that applies `resource_influx!` at `influx_time` intervals.
"""
function create_spatial_callbacks(t_inv, start_time, n_invaders, n_species, cutoff, n_comms, influx_amount, influx_time, n_complex, mucus, rng)
    invasion_cb = PeriodicCallback(
        integrator -> spatial_affect!(integrator, t_inv, start_time, n_invaders, n_species, cutoff, n_comms, rng), t_inv
    )
    if mucus
        n_influx_comms = n_comms -1
    else
        n_influx_comms = n_comms
    end
    resource_cb = PeriodicCallback(
        integrator -> resource_influx!(integrator, n_species, n_invaders, n_influx_comms, influx_amount, n_complex, rng), influx_time
    )
    return invasion_cb, resource_cb
end

