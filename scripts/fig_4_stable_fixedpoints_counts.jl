using Distributed
n_cores = length(Sys.cpu_info())
n_threaded_workers = floor(Int, n_cores / Base.Threads.nthreads())
# addprocs(n_threaded_workers - nworkers())
# @show nprocs()
using DrWatson
quickactivate("NeonateTriCorr")
using Pkg
Pkg.instantiate()

# @everywhere begin
using TravelingWaveSimulations, WilsonCowanModel, Simulation73
using Dates
using Makie
using GLMakie; ext_2d = "png"; GLMakie.activate!()
#using CairoMakie; ext_2d = "svg"; CairoMakie.activate!()
using AlgebraOfGraphics, ColorSchemes
using DataFrames
using ColorTypes
using NullclineAnalysis
using IterTools
using ThreadsX

include(srcdir("helpers.jl"))
include(srcdir("sweep.jl"))
include(srcdir("metrics.jl"))

# # loads failing_fp_arr and monotonic_fp_arr
# sub_A_sweep_lower_bound = 1.
# sub_A_sweep_upper_bound = 10.
# sub_A_range = sub_A_sweep_lower_bound..sub_A_sweep_upper_bound
include(scriptsdir("load/he_sweep_arrs.jl"))

if !@isdefined DEFAULT_NAME_MAPPING
    const DEFAULT_NAME_MAPPING = Dict("αE" => mods -> mods.α[1], "αI" => mods -> mods.α[2], "θE" => mods -> mods.θE, "a" => mods -> (@assert mods.aE == mods.firing_aI == mods.blocking_aI; mods.aE))
end

global_saved_len=15

nt_map(fn::Function, nt::NamedTuple{NS}) where NS = NamedTuple{NS}(map(fn, values(nt)))

function abbrev_count_label(x)
    if x >= 1000
        try
            "$(Int(x / 1000))K"
        catch
            "$(x / 1000)K"
        end
    else
        "$(Int(x))"
    end
end

function max_percentile_conditions(values; left_tail_pct=0.025, right_tail_pct=0.025)
    # Note: not symmetric tails, since values can be negative
    left_tail = maximum(values) * left_tail_pct
    right_tail = maximum(values) * (1 - right_tail_pct)
    middle_pct = 1 - (left_tail_pct + right_tail_pct)
    conditions = map(values) do value
        if value <= left_tail
            "min"
        elseif left_tail < value < right_tail
            "mid"
        elseif right_tail <= value
            "max"
        else
            throw(DomainError(value))
        end
    end
    return (conditions, left_tail, right_tail)
end

format_percent(x) = "$(round(Int, 100x))%"
format_tail(x) = "$(round(x, sigdigits=2))"

function subset_axes(axes_nt::NamedTuple{Names}, subset_range) where Names
    NamedTuple{Names}(
        [axis[axis .∈ subset_range] for axis in axes_nt]
    )
end

function plot_stable_fixedpoints_counts(mods; saved_lb, saved_ub, saved_len,
        subset_range=-Inf..Inf,
        name_mapping=DEFAULT_NAME_MAPPING,
        blocking_data = get_fp_arr_data("blocking", saved_lb, saved_ub, saved_len;
                mods=mods,
                name_mapping=name_mapping
            ), 
        monotonic_data = get_fp_arr_data("monotonic", saved_lb, saved_ub, saved_len;
                mods=mods,
                name_mapping=name_mapping
            ),
        blocking_fp_arr = blocking_data[1], 
        blocking_fp_axes = blocking_data[2], 
        blocking_prototype_name = blocking_data[4],
        monotonic_fp_arr = monotonic_data[1], 
        monotonic_fp_axes = monotonic_data[2], 
        monotonic_prototype_name = monotonic_data[4],
        fp_axes_nt = (
            FoI=subset_axes(blocking_fp_axes, subset_range),
            fire=subset_axes(monotonic_fp_axes, subset_range)
        ),
        fp_arr_nt = (
                FoI=blocking_fp_arr[
                    Aee=blocking_fp_axes.Aee .∈ subset_range,
                    Aei=blocking_fp_axes.Aei .∈ subset_range,
                    Aie=blocking_fp_axes.Aie .∈ subset_range,
                    Aii=blocking_fp_axes.Aii .∈ subset_range
                ],
                fire=monotonic_fp_arr[
                    Aee=monotonic_fp_axes.Aee .∈ subset_range,
                    Aei=monotonic_fp_axes.Aei .∈ subset_range,
                    Aie=monotonic_fp_axes.Aie .∈ subset_range,
                    Aii=monotonic_fp_axes.Aii .∈ subset_range
                ]
            ),
        nonl_types = keys(fp_axes_nt),
        prototype_name_nt = (
            fire=monotonic_prototype_name,
            FoI=blocking_prototype_name
        ),
        E_bounds = [0.05, 0.71],
        SI_bounds = [0.05, 0.71],
        arrows_step=0.05,
        # mods = (
        #     α=(0.4, 0.7), 
        #     firing_θI=0.2, blocking_θI=0.5, 
        #     n_lattice = 2,
        #     save_idxs=nothing, save_on=false, saveat=0.1
        # ),
        session_name = "fig_4_stable_fixedpoints_counts_nicolas",
        session_id = "$(Dates.now())",
        axis = (width = 800, height = 800),
        bar_theme = Theme(
            fontsize=56,
            strokewidth= 5.,
            Axis = (
                backgroundcolor = RGBA(1.,1.,1.,0.),
                leftspinevisible = true,
                rightspinevisible = false,
                bottomspinevisible = true,
                topspinevisible = false,
                xgridcolor = RGBA(1.,1.,1.,0.),
                ygridcolor = RGBA(1.,1.,1.,0.),
                strokewidth= 5.,
                ytickformat = xs -> abbrev_count_label.(xs)
            )
        ),
        nullcline_theme = Theme(
            fontsize=48,
            Axis = (
                backgroundcolor = RGBA(1.,1.,1.,0.),
                leftspinevisible = true,
                rightspinevisible = false,
                bottomspinevisible = true,
                topspinevisible = false,
                strokewidth=2.,
                xgridcolor = RGBA(1.,1.,1.,0.),
                ygridcolor = RGBA(1.,1.,1.,0.)
            ),
            Lines = (
                linewidth=4.0,
            ),
            Arrows = (
                arrowsize=10, lengthscale=0.017,
                linewidth=2,
                arrowcolor=:black, linecolor=:black,
                colormap=ColorSchemes.Greys_5,
                normalize=false
            ),
            Scatter = (
                markersize=27,
                strokewidth=1
            )
        ),
        plots_subdir = "$(session_name)_$(session_id)"
    )
    mkpath(plotsdir(plots_subdir))



    # # get only stable fixedpoints
    # stable_fp_arr_nt = NamedTuple{nonl_types}(
    #     map(nonl_types) do nonl_type
    #         filter_stable_fps(
    #             prototype_name_nt[nonl_type],
    #             mods,
    #             fp_arr_nt[nonl_type]
    #         )
    #     end
    # )

    # calculate seizure index of stable fixed points
    stablefp_SI_nt = NamedTuple{nonl_types}(
        map(nonl_types) do nonl_type
            prototype = get_prototype(prototype_name_nt[nonl_type])
            fp_arr = fp_arr_nt[nonl_type]
            fp_axes = fp_axes_nt[nonl_type]
            stablefp_SI = Array{Union{Float64,Missing}}(missing, size(fp_arr)..., 7)
            for (mx_idx, (nt, fps)) ∈ zip(CartesianIndices(fp_arr), enumerate_nda_dims(fp_arr, fp_axes))
                params = get_nullcline_params(prototype(; mods..., nt...))
                for fp_idx ∈ 1:length(fps)
                    fp = fps[fp_idx]
                    if fixedpoint_is_stable(fp, params)
                        stablefp_SI[mx_idx, fp_idx] = seizure_index(fp)
                    end
                end
            end
            stablefp_SI
        end
    )
    stablefp_E_nt = NamedTuple{nonl_types}(
        map(nonl_types) do nonl_type
            prototype = get_prototype(prototype_name_nt[nonl_type])
            fp_arr = fp_arr_nt[nonl_type]
            fp_axes = fp_axes_nt[nonl_type]
            stablefp_E = Array{Union{Float64,Missing}}(missing, size(fp_arr)..., 7)
            for (mx_idx, (nt, fps)) ∈ zip(CartesianIndices(fp_arr), enumerate_nda_dims(fp_arr, fp_axes))
                params = get_nullcline_params(prototype(; mods..., nt...))
                for fp_idx ∈ 1:length(fps)
                    fp = fps[fp_idx]
                    if fixedpoint_is_stable(fp, params)
                        stablefp_E[mx_idx, fp_idx] = first(fp)
                    end
                end
            end
            stablefp_E
        end
    )

    n_obs = length(stablefp_E_nt[nonl_types[1]])
    unrolled_E_df = DataFrame(
        nonl_type=repeat(nonl_types |> collect, inner=(n_obs,)),
        E = vcat([stablefp_E_nt[nonl][:] for nonl in nonl_types]...)
    )

    unrolled_SI_df = DataFrame(
        nonl_type=repeat(nonl_types |> collect, inner=(n_obs,)),
        SI = vcat([stablefp_SI_nt[nonl][:] for nonl in nonl_types]...)
    )
    dropmissing!(unrolled_E_df); dropmissing!(unrolled_SI_df)

    unrolled_SI_df.condition, SI_left, SI_right = max_percentile_conditions(unrolled_SI_df.SI)
    unrolled_E_df.condition, E_left, E_right = max_percentile_conditions(unrolled_E_df.E)

    condition_sorter = sorter("min", "mid", "max")

    with_theme(bar_theme) do

        SI_value_frequency = data(unrolled_SI_df) * histogram(bins=100) * mapping(:SI, layout=:nonl_type)

        fig = draw(SI_value_frequency; axis=merge(axis, (title="SI hist",)))
        save(plotsdir(plots_subdir, "stablefp_SI_histogram.$(ext_2d)"), fig)
    
        E_value_frequency = data(unrolled_E_df) * histogram(bins=100) * mapping(:E, layout=:nonl_type)
    
        fig = draw(E_value_frequency; axis=merge(axis, (title="E hist",)))
        save(plotsdir(plots_subdir, "stablefp_E_histogram.$(ext_2d)"), fig)

        SI_condition_frequency = data(unrolled_SI_df) * frequency() * mapping(:nonl_type)
        SI_stacks = SI_condition_frequency * mapping(color=:condition, stack=:condition => condition_sorter)
        fig = draw(SI_stacks; axis=merge(axis, (title="seizure index ($(format_tail(SI_left)), $(format_tail(SI_right)))",)))
        fig.figure.current_axis.x.xlabel = "inh. nonlinearity"
        tightlimits!.(filter(x->x isa Makie.Axis, contents(fig.figure.layout)))
        save(plotsdir(plots_subdir, "stablefp_SI_distribution.$(ext_2d)"), fig)

        ret_fig = fig

        E_condition_frequency = data(unrolled_E_df) * frequency() * mapping(:nonl_type)
        E_stacks = E_condition_frequency * mapping(color=:condition, stack=:condition => condition_sorter)
        fig = draw(E_stacks; axis=merge(axis, (title="excitatory activity ($(format_tail(E_left)), $(format_tail(E_right)))",)))
        fig.figure.current_axis.x.xlabel = "inh. nonlinearity"
        tightlimits!.(filter(x->x isa Makie.Axis, contents(fig.figure.layout)), Ref(Bottom()))
        save(plotsdir(plots_subdir, "stablefp_E_distribution.$(ext_2d)"), fig)

        ret_fig

    end # bar_theme

    # with_theme(nullcline_theme) do
    #     fig = Figure()
    #     fig[1,1] = ax = Makie.Axis(fig)
    #     xs = ones(Int, length(unrolled_SI_df.nonl_type))
    #     xs[unrolled_SI_df.nonl_type .== nonl_types[2]] .= 2
    #     violin!(ax, xs, unrolled_SI_df.SI, width=0.5)
    #     ax.xticks = [1,2]
    #     ax.xtickformat = xs -> ["monotonic", "failing"]
    #     fig
    # end


end # let

tau_name_mapping = merge(DEFAULT_NAME_MAPPING, Dict(:τE => m -> m.τ[1], :τI => m -> m.τ[2]))

# end # everywhere

# # Low A, unitary alpha
# let uniform_a = 50.,
#         mods=(
#             α=(1.0, 1.0), 
#             aE=uniform_a, firing_aI=uniform_a,
#             blocking_aI=uniform_a,
#             θE=0.125,
#             firing_θI=0.2, blocking_θI=0.5
#         ),
#         saved_lb=0.1, saved_ub=1.5, saved_len=10,
#         subset_range=saved_lb..saved_ub;
#     plot_stable_fixedpoints_counts(mods; saved_lb=saved_lb, saved_ub=saved_ub,
#         saved_len=saved_len, 
#         session_name="fig_4_stable_fixedpoints_counts_loA_unitAlpha"
#     )
# end

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
    saved_lb=1., saved_ub=20., saved_len=global_saved_len,
    subset_range=saved_lb..saved_ub;
    plot_stable_fixedpoints_counts(mods; saved_lb=saved_lb, saved_ub=saved_ub,
        saved_len=saved_len,
        name_mapping=tau_name_mapping,
        session_name="fig_4_stable_fixedpoints_counts_exhiA_unitAlpha",
        subset_range = subset_range
    )
end

# # Low-A non-unitary alpha
# let uniform_a = 50.,
#     mods=(
#         α=(0.4, 0.7), 
#         aE=uniform_a, firing_aI=uniform_a,
#         blocking_aI=uniform_a,
#         θE=0.125,
#         firing_θI=0.2, blocking_θI=0.5
#     ),
#     saved_lb=0.1, saved_ub=1.5, saved_len=global_saved_len,
#     subset_range=saved_lb..saved_ub;
# plot_stable_fixedpoints_counts(mods; saved_lb=saved_lb, saved_ub=saved_ub,
#     saved_len=saved_len, 
#     session_name="fig_4_stable_fixedpoints_counts_loA_nonunitAlpha"
# )
# end

# High-A non-unitary alpha
# let uniform_a = 5.,
#     mods=(
#         τ=(7.8, 7.8*4.4),
#         α=(0.4, 0.7), 
#         aE=uniform_a, firing_aI=uniform_a,
#         blocking_aI=uniform_a,
#         θE=1.25,
#         firing_θI=2.0, blocking_θI=5.0
#     ),
#     saved_lb=1., saved_ub=30., saved_len=global_saved_len,
#     subset_range=saved_lb..saved_ub;
# plot_stable_fixedpoints_counts(mods; saved_lb=saved_lb, saved_ub=saved_ub,
#     saved_len=saved_len, 
#     name_mapping=tau_name_mapping,
#     session_name="fig_4_stable_fixedpoints_counts_exhiA_nonunitAlpha",
#     subset_range = subset_range
# )
# end

