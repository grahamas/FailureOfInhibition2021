function get_fp_arr_data(nonl_name, lb, ub, len; 
        technical_mods=(
            n_lattice = 2, 
            save_idxs=nothing, 
            save_on=false, 
            saveat=0.1
        ),
        mods=(;),
        name_mapping=[],
        refresh_sweep_arrs=false
        )
    file, filename = produce_or_load(datadir(), 
        Dict("lb" => lb, 
             "ub"=> ub, 
             "len" => len, 
             [name => val_fn(mods) for (name, val_fn) in pairs(name_mapping)]...); 
        prefix = "$(nonl_name)_fp_arr",
        suffix = "bson",
        force = refresh_sweep_arrs,
        tag=false
    ) do c
        A_range = range(c["lb"], c["ub"], length=c["len"])
        nullcline_sweeping_mods = (Aee=A_range, Aei=A_range, Aie=A_range, Aii=A_range)
        nullcline_mods = merge(technical_mods, mods)
        prototype_name = "full_dynamics_$(nonl_name)"
        fp_arr = wcm_sweep_calculate_fixedpoints(
            prototype_name, 
            nullcline_mods,
            nullcline_sweeping_mods,
            ; 
            axis_length=100
        )
        return @dict(fp_arr, nullcline_mods, prototype_name)
    end
    @unpack fp_arr, nullcline_mods, prototype_name = file
    return (fp_arr, nullcline_mods, prototype_name)
end