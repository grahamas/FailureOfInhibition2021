using DrWatson
quickactivate("FailureOfInhibition2022")

# @everywhere begin
using TravelingWaveSimulations, WilsonCowanModel, Simulation73,
    NullclineAnalysis, Simulation73Plotting
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
include(scriptsdir("load/he_sweep_arrs.jl"))

function slice_at_dims_begin(A::NamedDimsArray, collapsed_dims::NTuple{N,Symbol}) where N
    dims_begins = axes.(Ref(A), collapsed_dims) .|> first
    dims_begin_index = NamedTuple{collapsed_dims |> Tuple}(dims_begins)
    getindex(A; dims_begin_index...)
end
function slice_at_dims_begin(A::NamedDimsArray{Names}, collapsed_dims::NTuple{N,Int}) where {Names,N}
    slice_at_dims_begin(A, Names[collapsed_dims |> collect] |> Tuple)
end

function slice_at_dims_end(A::NamedDimsArray, collapsed_dims::NTuple{N,Symbol}) where N
    dims_ends = axes.(Ref(A), collapsed_dims) .|> last
    dims_end_index = NamedTuple{collapsed_dims |> Tuple}(dims_ends)
    getindex(A; dims_end_index...)
end
function slice_at_dims_end(A::NamedDimsArray{Names}, collapsed_dims::NTuple{N,Int}) where {Names,N}
    slice_at_dims_end(A, Names[collapsed_dims |> collect] |> Tuple)
end

function count_stable_fps(fps::Vector{<:Vector}, prototype, mods, A_mods::NamedTuple)
    params = get_nullcline_params(prototype(; mods..., A_mods...))
    count(fixedpoint_is_stable.(fps, Ref(params)))
end
function count_stable_fps(saved_data::NamedTuple)
    prototype_name = saved_data.prototype_name
    prototype = get_prototype(prototype_name)
    fp_arr = saved_data.fp_arr
    fp_axes = saved_data.fp_axes
    mods = saved_data.nullcline_mods
    names = typeof(fp_arr).parameters[1]
    NamedDimsArray{names}(map(enumerate_nda_dims(fp_arr, fp_axes)) do (nt, fps)
        count_stable_fps(fps, prototype, mods, nt)
    end)
end

function explore_fp_counts(;
        sweep_lb, sweep_ub, sweep_len,  
        sweep_subset_range = -Inf..Inf,
        mods, name_mapping=DEFAULT_NAME_MAPPING,
        session_name = "explore_fp_counts",
        session_id = "$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))",
        plots_subdir = "$(session_name)_$(session_id)",
        figure_resolution=(800,800)
    )

    nonl_name_from_condition = Dict(
        Symbol("fire→fail") => "blocking",
        "fire" => "monotonic"
    )

    saved_data = Dict(
        condition => get_fp_arr_data(
            nonl_name, sweep_lb, sweep_ub, sweep_len;
            mods=mods,
            name_mapping=name_mapping
        ) for (condition, nonl_name) ∈ pairs(nonl_name_from_condition)
    )

    prototype_name_dct = Dict(
        condition => saved_data[condition].prototype_name for condition in keys(saved_data)
    )

    mkpath(plotsdir(plots_subdir))

    #fp_counts_by_condition = Dict(name => length.(getindex(saved_data, name).fp_arr) for name in keys(saved_data))
    stable_fp_counts_by_condition = Dict(
        name => count_stable_fps(saved_data[name]) for name in keys(saved_data)
    )

    figs_by_condition = map(keys(saved_data) |> collect) do condition
        fig = Figure(resolution=figure_resolution)
        fp_counts = stable_fp_counts_by_condition[condition]
        fp_axes = saved_data[condition].fp_axes
        layout = sweep_2D_slice_heatmaps(fig, fp_counts, fp_axes; 
            title="$(string(condition)) (end slices)", dimension_reduction=slice_at_dims_end
        )
        fig[1,1] = layout
        Makie.save(plotsdir(plots_subdir, "marginal_stable_fp_counts_$(condition).$(ext_2d)"), fig)
        fig
    end
    return figs_by_condition
end # function explore_fp_counts()

let uniform_a = 5.,
    sweep_params=(
        sweep_lb=1., sweep_ub=20., sweep_len=15,
        mods=(
            τ=(7.8, 7.8*4.4),
            α=(1.0, 1.0), 
            aE=uniform_a, firing_aI=uniform_a,
            blocking_aI=uniform_a,
            θE=1.25,
            firing_θI=2.0, blocking_θI=5.0
        )
    );
explore_fp_counts(; sweep_params...)
end