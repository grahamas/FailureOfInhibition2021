using DrWatson
@quickactivate "FailureOfInhibition2022"
using TravelingWaveSimulations, WilsonCowanModel, 
    Simulation73Plotting, Simulation73,
    NullclineAnalysis
using Dates
#using GLMakie; ext_2d = "png"; GLMakie.activate!()
using CairoMakie; ext_2d = "eps"; CairoMakie.activate!()

# loads blocking_fp_arr and monotonic_fp_arr
# include(scriptsdir("load/he_sweep_arrs.jl"))

include(srcdir("theme.jl"))
include(srcdir("sweep.jl"))
include(srcdir("plot_nullclines.jl"))
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")


function plot_nullclines_differing_balanced_state(;
        # sweep_lb, sweep_ub, sweep_len, 
        # sweep_subset_range = -Inf..Inf,
        # name_mapping=DEFAULT_NAME_MAPPING,
        mods, 
        arrows_step=0.05, 
        session_name = "fig_nullclines_differing_balanced_state",
        session_id = "$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))",
        plots_subdir = "$(session_name)_$(session_id)",
        figure_resolution=(1800,1000)
    )
    # saved_data = (
    #     blocking=get_fp_arr_data("blocking",sweep_lb,sweep_ub,sweep_len;
    #         mods=mods, name_mapping=name_mapping,
    #         subset_range=sweep_subset_range
    #     ),
    #     monotonic=get_fp_arr_data("monotonic",sweep_lb,sweep_ub,sweep_len;
    #         mods=mods, name_mapping=name_mapping,
    #         subset_range=sweep_subset_range
    #     )
    # )

    mkpath(plotsdir(plots_subdir))

    with_theme(nullcline_theme) do

    fig = Figure(resolution=figure_resolution)
    blocking_prototype = get_prototype("full_dynamics_blocking")

    # example_both_3sfp_idxs = findall(n_sfp_monotonic .== 3 .& n_sfp_blocking .== 3)
    # example_both_3sfp_idx = example_both_3sfp_idxs[end]
    # example_both_3sfp_coord = TravelingWaveSimulationsPlotting.get_coordinate(blocking_fp_axes, example_both_3sfp_idx)
    higher_aie = (Aee=19.0, Aei=13.0, Aie=24.0, Aii=6.0)
    higher_aie_nullcline_params = get_nullcline_params(blocking_prototype(; mods..., higher_aie...))
    @info "Higher Aie FPs: $(filter_stable_fps(calculate_fixedpoints(higher_aie_nullcline_params), higher_aie_nullcline_params))"
    fig[1,1] = plot_nullclines!(fig,
        get_nullcline_params(blocking_prototype(; mods..., higher_aie...));
        arrows_step=arrows_step,
        title="A[E → I] = $(higher_aie.Aie)"
    )

    mknt(; kwargs...) = NamedTuple{keys(kwargs)}(values(kwargs))

    lower_aie = merge(higher_aie, (Aie = 19.0,))
    lower_aie_nullcline_params = get_nullcline_params(blocking_prototype(; mods..., lower_aie...))
    @info "Lower Aie FPs: $(filter_stable_fps(calculate_fixedpoints(lower_aie_nullcline_params), lower_aie_nullcline_params))"
    fig[1,2] = plot_nullclines!(fig,
        lower_aie_nullcline_params;
        arrows_step=arrows_step,
        title="A[E → I] = $(lower_aie.Aie)"
    )
    hideydecorations!(fig[1,2] |> content)
    label_a = fig[1,1,TopLeft()] = Label(fig, "a", font=noto_sans_bold, halign=:left)
    label_b = fig[1,2,TopLeft()] = Label(fig, "b", font=noto_sans_bold, halign=:left)

    save(plotsdir(plots_subdir, "fig_nullclines_differing_balanced_state.$(ext_2d)"), fig)

    end # with_theme(nullcline_theme)
end # function


let uniform_a = 5.,
    sweep_params=(
        # sweep_lb=1., sweep_ub=20., sweep_len=15,
        mods=(
            τ=(7.8, 7.8*4.4),
            α=(1.0, 1.0), 
            aE=uniform_a, firing_aI=uniform_a,
            blocking_aI=uniform_a,
            θE=1.5,
            firing_θI=4.0, blocking_θI=8.0
        ),
    )
;
plot_nullclines_differing_balanced_state(;
    sweep_params...)
end
