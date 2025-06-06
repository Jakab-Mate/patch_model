function grab_ends(t, sim_length, transient_length, 
        gen, spec, counts, shannons, MCI)

    relaxation_length = sim_length - transient_length
    transient_end_indices = findall(t .<= transient_length .&& t .>= transient_length * 0.9)
    relaxation_end_indices = findall(t .<= sim_length .&& t .>= sim_length - relaxation_length * 0.1)

    gen_rel_end = mean(gen[relaxation_end_indices])
    spec_rel_end = mean(spec[relaxation_end_indices])
    counts_rel_end = mean(counts[relaxation_end_indices])
    shannons_rel_end = mean(shannons[relaxation_end_indices])
    MCI_rel_end = mean(MCI[relaxation_end_indices])

    gen_tr_end = mean(gen[transient_end_indices])
    spec_tr_end = mean(spec[transient_end_indices])
    counts_tr_end = mean(counts[transient_end_indices])
    shannons_tr_end = mean(shannons[transient_end_indices])
    MCI_tr_end = mean(MCI[transient_end_indices])

    gen_tr_rel = [gen_tr_end, gen_rel_end]
    spec_tr_rel = [spec_tr_end, spec_rel_end]
    counts_tr_rel = [counts_tr_end, counts_rel_end]
    shannons_tr_rel = [shannons_tr_end, shannons_rel_end]
    MCI_tr_rel = [MCI_tr_end, MCI_rel_end]

    return gen_tr_rel, spec_tr_rel, counts_tr_rel, shannons_tr_rel, MCI_tr_rel
end