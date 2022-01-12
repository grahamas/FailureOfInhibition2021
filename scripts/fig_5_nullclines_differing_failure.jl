using DrWatson
@quickactivate "FailureOfInhibition2021"
using TravelingWaveSimulations, WilsonCowanModel, 
    Simulation73Plotting, Simulation73,
    NullclineAnalysis
using Dates
#using GLMakie; ext_2d = "png"; GLMakie.activate!()
using CairoMakie; ext_2d = "svg"; CairoMakie.activate!()
using NamedDims, NamedDimsHelpers

include(srcdir("helpers.jl"))
include(srcdir("theme.jl"))
include(srcdir("sweep.jl"))
include(srcdir("metrics.jl"))
include(srcdir("plot_nullclines.jl"))
include(scriptsdir("load/he_sweep_arrs.jl"))

get_coordinate(names, axes, idx::CartesianIndex) = get_coordinate(names, axes, Tuple(idx))
function get_coordinate(names, axes, idx::Union{Tuple,Vector})
    @assert length(axes) == length(idx)
    NamedTuple{names}(map(zip(idx, axes)) do (i,ax)
        ax[i]
    end)
end

function plot_differing_failure(coord, mods; resolution, arrows_step)
    fig = Figure(resolution=resolution)
    fig[1,1] = plot_nullclines!(fig,
        get_nullcline_params(monotonic_prototype(; mods..., coord...));
        arrows_step=arrows_step,
        title = "firing"
    )
    fig[1,2] = plot_nullclines!(fig,
        get_nullcline_params(blocking_prototype(; mods..., both_3sfp_coord...));
        arrows_step=arrows_step,
        title = "firing → failing"
    )
    hideydecorations!(content(fig[1,2]))
    return fig
end

function fig_5_nullclines_differing_failure(; 
        mods, saved_lb, saved_ub, saved_len,
        session_name="fig_5_nullclines_differing_failure_lb=$(saved_lb)_ub=$(saved_ub)_len=$(saved_len)", 
        subset_range=saved_lb..saved_ub,
        name_mapping=DEFAULT_NAME_MAPPING,
        blocking_data = get_fp_arr_data("blocking", saved_lb, saved_ub, saved_len;
                mods=mods,
                name_mapping=name_mapping
            ), 
        monotonic_data = get_fp_arr_data("monotonic", saved_lb, saved_ub, saved_len;
                mods=mods,
                name_mapping=name_mapping
            ),
        blocking_fp_arr_unsubsetted::NDA = blocking_data[1], 
        blocking_fp_axes = blocking_data[2], 
        blocking_prototype_name = blocking_data[4],
        monotonic_fp_arr_unsubsetted::NDA = monotonic_data[1], 
        monotonic_fp_axes = monotonic_data[2], 
        monotonic_prototype_name = monotonic_data[4],
        fp_axes_nt = Dict(
            Symbol("fire→fail")=>subset_axes(blocking_fp_axes, subset_range),
            :fire=>subset_axes(monotonic_fp_axes, subset_range)
        ),
        fp_arr_nt = Dict(
                Symbol("fire→fail")=>blocking_fp_arr_unsubsetted[
                    Aee=blocking_fp_axes.Aee .∈ subset_range,
                    Aei=blocking_fp_axes.Aei .∈ subset_range,
                    Aie=blocking_fp_axes.Aie .∈ subset_range,
                    Aii=blocking_fp_axes.Aii .∈ subset_range
                ],
                :fire=>monotonic_fp_arr_unsubsetted[
                    Aee=monotonic_fp_axes.Aee .∈ subset_range,
                    Aei=monotonic_fp_axes.Aei .∈ subset_range,
                    Aie=monotonic_fp_axes.Aie .∈ subset_range,
                    Aii=monotonic_fp_axes.Aii .∈ subset_range
                ]
            ),
        nonl_types = keys(fp_axes_nt) |> Tuple,
        prototype_name_nt = Dict(
            :fire=>monotonic_prototype_name,
            Symbol("fire→fail")=>blocking_prototype_name
        ),
        arrows_step=0.05,
        session_id = "$(Dates.now())",
        axis = (width = 1800, height = 1000),
        figure_resolution = (1800,1000),
        plots_subdir = "$(session_name)_$(session_id)"
    ) where {Names,NDA<:NamedDimsArray{Names}}
    mkpath(plotsdir(plots_subdir))

    blocking_fp_arr = fp_arr_nt[Symbol("fire→fail")];
    monotonic_fp_arr = fp_arr_nt[:fire];

    blocking_fp_axes = fp_axes_nt[Symbol("fire→fail")];
    monotonic_fp_axes = fp_axes_nt[:fire];

    n_sfp_blocking = count_stable_fps(blocking_fp_arr, blocking_fp_axes, blocking_prototype, mods);
    n_sfp_monotonic = count_stable_fps(monotonic_fp_arr, monotonic_fp_axes, monotonic_prototype, mods);

    @show extrema(n_sfp_blocking)
    @show extrema(n_sfp_monotonic)

    target_idxs = findall(n_sfp_blocking .== 3);
    isempty(target_idxs) && error("No sim with blocking 3 SFP found.")
    target_idx = rand(target_idxs)
    target_coord = get_coordinate(Names, blocking_fp_axes, target_idx)
    @show "Both: $target_coord"

    with_theme(nullcline_theme) do
        differing_failure_fig = plot_differing_failure(target_coord, mods; resolution=figure_resolution, arrows_step=arrows_step)

        save(plotsdir(plots_subdir, "fig_5_nullclines_differing_failure.$(ext_2d)"), differing_failure_fig)
    end # with_theme(nullcline_theme)
end # function fig_5


# High-A unitary alpha
let uniform_a = 5.,
    mods=(
        τ=(7.8, 7.8*4.4),
        α=(1.0, 1.0), 
        aE=uniform_a, firing_aI=uniform_a,
        blocking_aI=uniform_a,
        θE=1.25,
        firing_θI=2.0, blocking_θI=5.0
    ),
    saved_lb=1., saved_ub=100., saved_len=10,
    subset_range=saved_lb..saved_ub;
    fig_5_nullclines_differing_failure(; mods=mods, 
        saved_lb=saved_lb, saved_ub=saved_ub, saved_len=saved_len
    )
end