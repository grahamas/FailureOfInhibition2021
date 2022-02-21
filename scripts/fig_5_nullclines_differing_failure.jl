using DrWatson
@quickactivate "FailureOfInhibition2022"
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

function plot_differing_failure(coord, mods, prototypes, titles; resolution, arrows_step)
    fig = Figure(resolution=resolution)
    for (i,(key,title)) in enumerate(titles)
        fig[1,i] = plot_nullclines!(fig,
            get_nullcline_params(prototypes[key](; mods..., coord...));
            arrows_step=arrows_step,
            title = title
        )
    end
    hideydecorations!.([content(fig[1,i]) for i ∈ 2:length(prototypes)])
    return fig
end

function ⋉(needles, haystack)
    [needle ∈ haystack for needle in needles]
end

function fig_5_nullclines_differing_failure(; 
        mods, saved_lb, saved_ub, saved_len,
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
        fp_axes_dct = Dict(
            Symbol("fire→fail")=>subset_axes(blocking_fp_axes, subset_range),
            :fire=>subset_axes(monotonic_fp_axes, subset_range)
        ),
        fp_arr_dct = Dict(
                Symbol("fire→fail")=>blocking_fp_arr_unsubsetted[
                    Aee=blocking_fp_axes.Aee ⋉ subset_range,
                    Aei=blocking_fp_axes.Aei ⋉ subset_range,
                    Aie=blocking_fp_axes.Aie ⋉ subset_range,
                    Aii=blocking_fp_axes.Aii ⋉ subset_range
                ],
                :fire=>monotonic_fp_arr_unsubsetted[
                    Aee=monotonic_fp_axes.Aee ⋉ subset_range,
                    Aei=monotonic_fp_axes.Aei ⋉ subset_range,
                    Aie=monotonic_fp_axes.Aie ⋉ subset_range,
                    Aii=monotonic_fp_axes.Aii ⋉ subset_range
                ]
            ),
        nonl_types = keys(fp_axes_dct) |> Tuple,
        prototype_name_dct = Dict(
            :fire=>monotonic_prototype_name,
            Symbol("fire→fail")=>blocking_prototype_name
        ),
        title_pairs = [
            :fire=>"fire",
            Symbol("fire→fail")=>"firing → failing"
        ],
        arrows_step=0.05,
        axis = (width = 1800, height = 1000),
        figure_resolution = (1800,1000),
        session_name="fig_5_nullclines_differing_failure_lb=$(saved_lb)_ub=$(saved_ub)_len=$(saved_len)", 
        session_id = "$(Dates.now())",
        plots_subdir = "$(session_name)_$(session_id)"
    ) where {Names,NDA<:NamedDimsArray{Names}}
    if !isdir(plotsdir(plots_subdir))
        mkpath(plotsdir(plots_subdir))
    end

    prototype_dct = Dict(
        sym => get_prototype(prototype_name_dct[sym]) for sym in keys(fp_axes_dct)
    )
    n_sfp = Dict(
        sym => count_stable_fps(fp_arr_dct[sym], fp_axes_dct[sym], prototype_dct[sym], mods) for sym ∈ keys(fp_axes_dct)
    )
    n_fps = Dict(
        sym => length.(fp_arr_dct[sym]) for sym ∈ keys(fp_axes_dct)
    )
    n_relevant_fps = n_fps
    conds = Dict(
        Symbol("fire→fail") => n -> n >= 7,
        Symbol("fire") => n -> n >= 5 
    )
    conds_strs = Dict(
        Symbol("fire→fail") => ">= 7",
        Symbol("fire") => ">= 5"
    )

    for name in keys(n_relevant_fps)
        @info "$(name) extrema: $(extrema(n_relevant_fps[name]))"
    end

    conds_satisfied_dct = Dict(
        sym => conds[sym].(n_relevant_fps[sym]) for sym in keys(conds)
    )
    and(x,y) = x && y
    target_idxs = findall(reduce((x,y) -> and.(x,y), values(conds_satisfied_dct)));
    isempty(target_idxs) && error("No sim with $(["$sym $conds" for (sym,conds) ∈ pairs(conds_strs)]) relevant FPs found.")
    target_idx = rand(target_idxs)
    target_coord = get_coordinate(Names, blocking_fp_axes, target_idx)
    @show "Both: $target_coord"

    with_theme(nullcline_theme) do
        differing_failure_fig = plot_differing_failure(target_coord, mods, prototype_dct, title_pairs; resolution=figure_resolution, arrows_step=arrows_step)

        save(plotsdir(plots_subdir, "fig_5_nullclines_differing_failure_fps_5_7.$(ext_2d)"), differing_failure_fig)
    end # with_theme(nullcline_theme)
end # function fig_5


# High-A unitary alpha
let uniform_a = 5.,
    mods=(
        τ=(7.8, 7.8),
        α=(1.0,1.0), 
        aE=uniform_a, firing_aI=uniform_a,
        blocking_aI=uniform_a,
        θE=1.5,
        firing_θI=4.0, blocking_θI=8.0
    ),
    saved_lb=1., saved_ub=20., saved_len=15,
    subset_range=saved_lb..saved_ub,
    name_mapping = DEFAULT_NAME_MAPPING,
    session_name="fig_5_nullclines_differing_failure_fps_5_7_lb=$(saved_lb)_ub=$(saved_ub)_len=$(saved_len)", 
    session_id = "$(Dates.now())",
    plots_subdir = "$(session_name)_$(session_id)";
    if !isdir(plotsdir(plots_subdir))
        mkpath(plotsdir(plots_subdir))
    end
    fig_5_nullclines_differing_failure(; mods=mods, 
        saved_lb=saved_lb, saved_ub=saved_ub, saved_len=saved_len,
        session_name=session_name, session_id=session_id, plots_subdir=plots_subdir
    )
end