
if !@isdefined DEFAULT_NAME_MAPPING
    const DEFAULT_NAME_MAPPING = Dict("αE" => mods -> mods.α[1], "αI" => mods -> mods.α[2], "θE" => mods -> mods.θE, "a" => mods -> (@assert mods.aE == mods.firing_aI == mods.blocking_aI; mods.aE), :τE => m -> m.τ[1], :τI => m -> m.τ[2])
end

function get_fp_arr_data(nonl_name, lb, ub, len; 
        technical_mods=(
            n_lattice = 2, 
            save_idxs=nothing, 
            save_on=false, 
            saveat=0.1
        ),
        mods=(;),
        name_mapping=[],
        refresh_sweep_arrs=false,
        subset_range=-Inf..Inf
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
        fp_arr = threaded_sweep_calculate_fixedpoints(
            prototype_name, 
            nullcline_mods,
            nullcline_sweeping_mods,
            ; 
            axis_length=100
        )
        fp_axes = nullcline_sweeping_mods
        return @dict(fp_arr, fp_axes, nullcline_mods, prototype_name)
    end
    @unpack fp_arr, fp_axes, nullcline_mods, prototype_name = file
    fp_arr, fp_axes = subset_nda_dims(fp_arr, fp_axes, subset_range)
    return (fp_arr=fp_arr, fp_axes=fp_axes, nullcline_mods=nullcline_mods, prototype_name=prototype_name)
end