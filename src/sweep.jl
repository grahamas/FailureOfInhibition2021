# using DrWatson
# using TravelingWaveSimulations, WilsonCowanModel,
#     Simulation73Plotting, Simulation73
# using Dates
# using Makie
# using AxisIndices, IterTools

using ProgressMeter

function sweep_calculate_fixedpoints(prototype_name::String, static_mods, dynamic_mods::NamedTuple{NAMES}; axis_length::Integer=150) where {NAMES}
    if :n_lattice ∉ NAMES
        # @warn "Don't forget to make n_lattice small!"
    end
    prototype = get_prototype(prototype_name)
    sweep_axes = values(dynamic_mods)

    us = range(0., 1., length=axis_length)
    vs = copy(us)
    dus = Array{Float64,2}(undef, length(us), length(vs))

    sweep = product(sweep_axes...) |> collect
    progress = Progress(length(sweep); dt=1, desc="Calculating fixed points array...", showspeed=true)

    # FIXME potential parallelization target -- but arrays variable len
    NamedDimsArray{NAMES}(map(sweep) do (sweeping_vals...)
        sweeping_mods = NamedTuple{NAMES}(sweeping_vals...)
        sim = prototype(; static_mods..., sweeping_mods...)
        params = get_nullcline_params(sim.model)
        fps = calculate_fixedpoints!(dus, [us, vs], (WilsonCowanModel.wcm_du_defn, WilsonCowanModel.wcm_dv_defn), params, ((0.,1.), (0.,1.)))
        if length(fps) > 7
            @warn "$(sweeping_vals) have $(length(fps)) fixed points."
        end
        ProgressMeter.next!(progress)
        fps
    end)
end

function distributed_sweep_calculate_fixedpoints(prototype_name::String, static_mods, dynamic_mods::NamedTuple{NAMES}; axis_length::Integer=150) where {NAMES}
    if :n_lattice ∉ NAMES
        # @warn "Don't forget to make n_lattice small!"
    end
    prototype = get_prototype(prototype_name)
    sweep_axes = values(dynamic_mods)

    sweep = product(sweep_axes...) |> collect
    progress = Progress(length(sweep); dt=1, desc="Calculating fixed points array...", showspeed=true)

    progress_lock = ReentrantLock()
    prototype_lock = ReentrantLock()

    # FIXME potential parallelization target -- but arrays variable len
    NamedDimsArray{NAMES}(pmap(sweep) do (sweeping_vals...)
        sweeping_mods = NamedTuple{NAMES}(sweeping_vals...)

        sim, ax_len = lock(prototype_lock) do
            prototype(; static_mods..., sweeping_mods...), axis_length
        end

        us = range(0., 1., length=axis_length)
        vs = copy(us)
        dus = Array{Float64,2}(undef, length(us), length(vs))

        params = get_nullcline_params(sim.model)
        fps = calculate_fixedpoints!(dus, [us, vs], (WilsonCowanModel.wcm_du_defn, WilsonCowanModel.wcm_dv_defn), params, ((0.,1.), (0.,1.)))
        if length(fps) > 7
            @warn "$(sweeping_vals) have $(length(fps)) fixed points."
        end
        lock(progress_lock) do
            ProgressMeter.next!(progress)
        end
        fps
    end)
end

function threaded_sweep_calculate_fixedpoints(prototype_name::String, static_mods, dynamic_mods::NamedTuple{NAMES}; axis_length::Integer=150) where {NAMES}
    if :n_lattice ∉ NAMES
        # @warn "Don't forget to make n_lattice small!"
    end
    prototype = get_prototype(prototype_name)
    sweep_axes = values(dynamic_mods)

    sweep = product(sweep_axes...) |> collect
    progress = Progress(length(sweep); dt=1, desc="Calculating fixed points array...", showspeed=true)

    progress_lock = ReentrantLock()
    prototype_lock = ReentrantLock()

    # FIXME potential parallelization target -- but arrays variable len
    NamedDimsArray{NAMES}(ThreadsX.map(sweep) do (sweeping_vals...)
        sim, ax_len, fns = lock(prototype_lock) do
            sweeping_mods = NamedTuple{NAMES}(sweeping_vals...)
            (prototype(; static_mods..., sweeping_mods...), axis_length, (WilsonCowanModel.wcm_du_defn, WilsonCowanModel.wcm_dv_defn))
        end

        us = range(0., 1., length=ax_len)
        vs = copy(us)
        dus = Array{Float64,2}(undef, length(us), length(vs))

        params = get_nullcline_params(sim.model)
        fps = calculate_fixedpoints!(dus, [us, vs], fns, params, ((0.,1.), (0.,1.)))
        if length(fps) > 7
            @warn "$(sweeping_vals) have $(length(fps)) fixed points."
        end
        lock(progress_lock) do
            ProgressMeter.next!(progress)
        end
        fps
    end)
end