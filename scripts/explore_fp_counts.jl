using DrWatson
quickactivate("FailureOfInhibition2021")

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

function explore_fp_counts(;
        sweep_lb, sweep_ub, sweep_len,  
        sweep_subset_range = -Inf..Inf,
        mods, name_mapping=DEFAULT_NAME_MAPPING,
        session_name = "explore_fp_counts",
        session_id = "$(Dates.now())",
        plots_subdir = "$(session_name)_$(session_id)",
        figure_resolution=(800,800)
    )

    saved_data = (
        fire_fail = get_fp_arr_data(
            "blocking", sweep_lb, sweep_ub, sweep_len;
            mods=mods,
            name_mapping=name_mapping
        ),
        fire = get_fp_arr_data(
            "monotonic", sweep_lb, sweep_ub, sweep_len;
            mods=mods,
            name_mapping=name_mapping
        )
    )

    mkpath(plotsdir(plots_subdir))

    fp_counts_by_condition = Dict(name => length.(getproperty(saved_data, name).fp_arr) for name in keys(saved_data))

    figs_by_condition = map(keys(saved_data)) do condition
        fig = Figure(resolution=figure_resolution)
        fp_counts = fp_counts_by_condition[condition]
        fp_axes = saved_data[condition].fp_axes
        layout = sweep_2D_slice_heatmaps(fig, fp_counts, fp_axes; 
            title=string(condition)
        )
        fig[1,1] = layout
        Makie.save(plotsdir(plots_subdir, "marginal_fp_counts_$(condition).$(ext_2d)"), fig)
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