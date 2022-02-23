using DrWatson
@quickactivate "FailureOfInhibition2022"
using TravelingWaveSimulations, WilsonCowanModel, 
    Simulation73Plotting, Simulation73,
    NullclineAnalysis
using Dates
#using GLMakie; ext_2d = "png"; GLMakie.activate!()
using CairoMakie; ext_2d = "png"; CairoMakie.activate!()
using NamedDims, NamedDimsHelpers

include(srcdir("helpers.jl"))
include(srcdir("theme.jl"))
include(srcdir("sweep.jl"))
include(srcdir("metrics.jl"))
include(srcdir("plot_nullclines.jl"))
include(scriptsdir("load/he_sweep_arrs.jl"))

noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")

tuple_str(tup) = join(["$key=$val" for (key,val) in pairs(tup)], "_")

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
    label_a = fig[1,1,TopLeft()] = Label(fig, "a", font=noto_sans_bold, halign=:left)
    label_b = fig[1,2,TopLeft()] = Label(fig, "b", font=noto_sans_bold, halign=:left)
    return fig
end

function ⋉(needles, haystack)
    [needle ∈ haystack for needle in needles]
end

# function fig_nullclines_differing_failure(; 
#         mods, saved_lb, saved_ub, saved_len,
#         subset_range=saved_lb..saved_ub,
#         name_mapping=DEFAULT_NAME_MAPPING,
#         blocking_data = get_fp_arr_data("blocking", saved_lb, saved_ub, saved_len;
#                 mods=mods,
#                 name_mapping=name_mapping
#             ), 
#         monotonic_data = get_fp_arr_data("monotonic", saved_lb, saved_ub, saved_len;
#                 mods=mods,
#                 name_mapping=name_mapping
#             ),
#         blocking_fp_arr_unsubsetted::NDA = blocking_data[1], 
#         blocking_fp_axes = blocking_data[2], 
#         blocking_prototype_name = blocking_data[4],
#         monotonic_fp_arr_unsubsetted::NDA = monotonic_data[1], 
#         monotonic_fp_axes = monotonic_data[2], 
#         monotonic_prototype_name = monotonic_data[4],
#         fp_axes_dct = Dict(
#             Symbol("fire→fail")=>subset_axes(blocking_fp_axes, subset_range),
#             :fire=>subset_axes(monotonic_fp_axes, subset_range)
#         ),
#         fp_arr_dct = Dict(
#                 Symbol("fire→fail")=>blocking_fp_arr_unsubsetted[
#                     Aee=blocking_fp_axes.Aee ⋉ subset_range,
#                     Aei=blocking_fp_axes.Aei ⋉ subset_range,
#                     Aie=blocking_fp_axes.Aie ⋉ subset_range,
#                     Aii=blocking_fp_axes.Aii ⋉ subset_range
#                 ],
#                 :fire=>monotonic_fp_arr_unsubsetted[
#                     Aee=monotonic_fp_axes.Aee ⋉ subset_range,
#                     Aei=monotonic_fp_axes.Aei ⋉ subset_range,
#                     Aie=monotonic_fp_axes.Aie ⋉ subset_range,
#                     Aii=monotonic_fp_axes.Aii ⋉ subset_range
#                 ]
#             ),
#         nonl_types = keys(fp_axes_dct) |> Tuple,
#         prototype_name_dct = Dict(
#             :fire=>monotonic_prototype_name,
#             Symbol("fire→fail")=>blocking_prototype_name
#         ),
#         title_pairs = [
#             :fire=>"fire",
#             Symbol("fire→fail")=>"firing → failing"
#         ],
#         arrows_step=0.05,
#         axis = (width = 1800, height = 1000),
#         figure_resolution = (1800,1000),
#         session_name="fig_nullclines_differing_failure_lb=$(saved_lb)_ub=$(saved_ub)_len=$(saved_len)", 
#         session_id = "$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))",
#         plots_subdir = "$(session_name)_$(session_id)"
#     ) where {Names,NDA<:NamedDimsArray{Names}}
#     if !isdir(plotsdir(plots_subdir))
#         mkpath(plotsdir(plots_subdir))
#     end

#     prototype_dct = Dict(
#         sym => get_prototype(prototype_name_dct[sym]) for sym in keys(fp_axes_dct)
#     )
#     n_sfp = Dict(
#         sym => count_stable_fps(fp_arr_dct[sym], fp_axes_dct[sym], prototype_dct[sym], mods) for sym ∈ keys(fp_axes_dct)
#     )
#     n_sfp_target = Dict(
#         Symbol("fire→fail") => 3
#     )

#     for name in keys(n_sfp)
#         @info "$(name) extrema: $(extrema(n_sfp[name]))"
#     end

#     equals_target_dct = Dict(
#         sym => n_sfp[sym] .== n_sfp_target[sym] for sym in keys(n_sfp_target)
#     )
#     and(x,y) = x && y
#     target_idxs = findall(reduce((x,y) -> and.(x,y), values(equals_target_dct)));
#     isempty(target_idxs) && error("No sim with $(["$sym == $target" for (sym,target) ∈ pairs(n_sfp_target)]) SFP found.")
#     target_idx = rand(target_idxs)
#     target_coord = get_coordinate(Names, blocking_fp_axes, target_idx)
#     @show "Both: $target_coord"

#     with_theme(nullcline_theme) do
#         differing_failure_fig = plot_differing_failure(target_coord, mods, prototype_dct, title_pairs; resolution=figure_resolution, arrows_step=arrows_step)

#         save(plotsdir(plots_subdir, "fig_nullclines_differing_failure.$(ext_2d)"), differing_failure_fig)
#     end # with_theme(nullcline_theme)
# end # function fig


# High-A unitary alpha
let uniform_a = 5.,
    mods=(
        τ=(7.8, 7.8*4.4),
        α=(0.4,0.7), 
        aE=uniform_a, firing_aI=uniform_a,
        blocking_aI=uniform_a,
        θE=1.5,
        firing_θI=4.0, blocking_θI=8.0
    ),
    saved_lb=1., saved_ub=20., saved_len=15,
    subset_range=saved_lb..saved_ub,
    name_mapping = DEFAULT_NAME_MAPPING,
    session_name="fig_nullclines_differing_failure_lb=$(saved_lb)_ub=$(saved_ub)_len=$(saved_len)", 
    session_id = "$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))",
    plots_subdir = "$(session_name)_$(session_id)";
    if !isdir(plotsdir(plots_subdir))
        mkpath(plotsdir(plots_subdir))
    end
    # fig_nullclines_differing_failure(; mods=mods, 
    #     saved_lb=saved_lb, saved_ub=saved_ub, saved_len=saved_len,
    #     session_name=session_name, session_id=session_id, plots_subdir=plots_subdir
    # )
    prototype_name_dct = Dict(
            :fire=>"full_dynamics_monotonic",
            Symbol("fire→fail")=>"full_dynamics_blocking"
        )
    prototype_dct = Dict(
        sym => get_prototype(name) for (sym, name) in pairs(prototype_name_dct)
    )
    title_pairs = [
            :fire=>"fire",
            Symbol("fire→fail")=>"firing → failing"
    ]
    with_theme(nullcline_theme) do 
        mod = (
            Aee=10., Aei=7., Aie=20., Aii=6., θE=1.5,
            firing_θI=4.0, blocking_θI=12.0, τ=(7.8, 7.8),
            α=(.4,.7)
        )
        mod1 = (Aee=11.0,)
        mod2 = (Aee=12.0,)
        mod3 = (Aee=13.0,)
        mod4 = (Aee=11.5,)

        @show mod
        mod_fig = plot_differing_failure(merge(mod, mod), mods, prototype_dct, title_pairs; resolution=(1800,1000), arrows_step=0.05)
        save(plotsdir(plots_subdir, "fig_nullclines_differing_failure_$(tuple_str(mod)).$(ext_2d)"), mod_fig)

        mod1_fig = plot_differing_failure(merge(mod, mod1), mods, prototype_dct, title_pairs; resolution=(1800,1000), arrows_step=0.05)
        save(plotsdir(plots_subdir, "fig_nullclines_differing_failure_$(tuple_str(mod1)).$(ext_2d)"), mod1_fig)

        mod2_fig = plot_differing_failure(merge(mod, mod2), mods, prototype_dct, title_pairs; resolution=(1800,1000), arrows_step=0.05)
        save(plotsdir(plots_subdir, "fig_nullclines_differing_failure_$(tuple_str(mod2)).$(ext_2d)"), mod2_fig)

        mod3_fig = plot_differing_failure(merge(mod, mod3), mods, prototype_dct, title_pairs; resolution=(1800,1000), arrows_step=0.05)
        save(plotsdir(plots_subdir, "fig_nullclines_differing_failure_$(tuple_str(mod3)).$(ext_2d)"), mod3_fig)

        mod4_fig = plot_differing_failure(merge(mod, mod4), mods, prototype_dct, title_pairs; resolution=(1800,1000), arrows_step=0.05)
        save(plotsdir(plots_subdir, "fig_nullclines_differing_failure_$(tuple_str(mod4)).$(ext_2d)"), mod4_fig)
    end
end