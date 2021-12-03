using DrWatson
@quickactivate "FailureOfInhibition2021"
using TravelingWaveSimulations, WilsonCowanModel, 
    TravelingWaveSimulationsPlotting, 
    Simulation73Plotting, Simulation73
using Dates
#using GLMakie; ext_2d = "png"; GLMakie.activate!()
using CairoMakie; ext_2d = "svg"; CairoMakie.activate!()
using AxisIndices 

# loads blocking_fp_arr and monotonic_fp_arr
sub_A_sweep_lower_bound = 0.5
sub_A_sweep_upper_bound = 1.5
sub_A_range = sub_A_sweep_lower_bound..sub_A_sweep_upper_bound
include(scriptsdir("load/he_sweep_arrs.jl"))

include(srcdir("theme.jl"))


let blocking_fp_arr = blocking_fp_arr[Aee=sub_A_range,Aei=sub_A_range,Aie=sub_A_range,Aii=sub_A_range], 
    monotonic_fp_arr = monotonic_fp_arr[Aee=sub_A_range,Aei=sub_A_range,Aie=sub_A_range,Aii=sub_A_range],
    blocking_fp_count_arr = length.(blocking_fp_arr), 
    monotonic_fp_count_arr = length.(monotonic_fp_arr),
    blocking_prototype_name = "full_dynamics_blocking",
    monotonic_prototype_name = "full_dynamics_monotonic",
    arrows_step=0.05,
    smods = (
        α=(0.4, 0.7), 
        firing_θI=0.2, blocking_θI=0.5, 
        n_lattice = 2,
        save_idxs=nothing, save_on=false, saveat=0.1
    ),
    session_name = "fig_5_nullclines_differing_failure",
    session_id = "$(Dates.now())",
    plots_subdir = "$(session_name)_$(session_id)",
    figure_resolution=(1800,1000)
;
mkpath(plotsdir(plots_subdir))

monotonic_prototype = get_prototype(monotonic_prototype_name)
blocking_prototype = get_prototype(blocking_prototype_name)

blocking_stable_fp_arr = filter_stable_fps(blocking_prototype_name, smods, blocking_fp_arr)
monotonic_stable_fp_arr = filter_stable_fps(monotonic_prototype_name, smods, monotonic_fp_arr)

n_sfp_blocking = count_stable_fps(blocking_prototype_name, smods, blocking_fp_arr)
n_sfp_monotonic = count_stable_fps(monotonic_prototype_name, smods, monotonic_fp_arr)

with_theme(nullcline_theme) do

both_3sfp_idxs = findall(n_sfp_monotonic .== 3 .& n_sfp_blocking .== 3)
both_3sfp_idx = both_3sfp_idxs[end]
both_3sfp_coord = TravelingWaveSimulationsPlotting.get_coordinate(blocking_fp_arr, both_3sfp_idx)
@show "Both: $both_3sfp_coord"

both_3sfp_fig = Figure(resolution=figure_resolution)
both_3sfp_fig[1,1] = plot_nullclines!(both_3sfp_fig,
    get_nullcline_params(monotonic_prototype(; smods..., both_3sfp_coord...));
    arrows_step=arrows_step,
    title = "firing"
)
both_3sfp_fig[1,2] = plot_nullclines!(both_3sfp_fig,
    get_nullcline_params(blocking_prototype(; smods..., both_3sfp_coord...));
    arrows_step=arrows_step,
    title = "firing → failing"
)
hideydecorations!(content(both_3sfp_fig[1,2]))


TravelingWaveSimulationsPlotting._save!(both_3sfp_fig, "fig_5_nullclines_differing_failure.$(ext_2d)"; plots_subdir=plots_subdir)

end # with_theme(nullcline_theme)
end # let

