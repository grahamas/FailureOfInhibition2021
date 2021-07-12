using Base: NamedTuple
using DrWatson
@quickactivate "FailureOfInhibition2021"

# using NeuralModels: AxisIndices

using Base.Iterators
using Optim, Statistics, AxisIndices

using CairoMakie,Makie

using NeuralModels

using Memoize

include(srcdir("axisarray.jl"))
include(srcdir("estimate_params.jl"))
include(srcdir("theme.jl"))

if !@isdefined(force_results_calc) || force_results_calc
    μ_on=0.; σ_on=0.04:0.02:0.24; 
    μ_off=0.5:0.05:1.; σ_off=0.04:0.02:0.24    
    results = fit_range_of_thresholds(;
        μ_on=μ_on, σ_on=σ_on,
        μ_off=μ_off, σ_off=σ_off
    )
    force_results_calc = false
end

force_results_calc = let results = results;

results_df = results_axisarray_to_df(results; model_names=[:my, :meijer])#, :kim])

model_string_names = Dict(
    :my => "DoS",
    :meijer => "Gauss",
    :kim => "Kim"
)

with_theme(nullcline_theme) do
    rmse_plt = data(results_df) * mapping(
        :samples => (s -> s.off.μ - s.on.μ) => "Δμ",
        :fits => Optim.minimum => "RMSE",
        color=(:model => m -> model_string_names[m])
    ) * (visual(Scatter; markersize=5, strokewidth=0)+ linear(linewidth=3))
    wide_idσ_result = results[σ_on=σ_on[begin],
        μ_off=1., σ_off=σ_off[begin]]
    wide_transσ_result = results[σ_on=σ_on[3end÷4],
        μ_off=1., σ_off=σ_off[begin]]
    narrow_transσ_result = results[σ_on=σ_on[3end÷4],
        μ_off=0.5, σ_off=σ_off[begin]]
    drawn_fig = draw(rmse_plt)#; axis)
    fig = drawn_fig.figure
    fig.scene.resolution[] = (1100,1000)
    #axislegend(Makie.Axis(first(drawn_fig)); position=:lt)

    wide_idσ_Δμ = wide_idσ_result.samples.off.μ - wide_idσ_result.samples.on.μ
    wide_idσ_title = "Δμ = $(round(wide_idσ_Δμ,digits=2)), σₒₙ=$(round(wide_idσ_result.samples.on.σ, digits=2))"
    fig[0,1] = wide_idσ_ax = Makie.Axis(fig, title=wide_idσ_title)
    @info """$wide_idσ_title
        DoS RMSE: $(Optim.minimum(wide_idσ_result.fits.my))
        Gaussian RMSE: $(Optim.minimum(wide_idσ_result.fits.meijer))"""

    wide_transσ_Δμ = wide_transσ_result.samples.off.μ - wide_transσ_result.samples.on.μ
    wide_transσ_title = "Δμ = $(round(wide_transσ_Δμ,digits=2)), σₒₙ=$(round(wide_transσ_result.samples.on.σ, digits=2))"
    fig[1,0] = wide_transσ_ax = Makie.Axis(fig, title=wide_transσ_title)
    @info """$wide_transσ_title
        DoS RMSE: $(Optim.minimum(wide_transσ_result.fits.my))
        Gaussian RMSE: $(Optim.minimum(wide_transσ_result.fits.meijer))"""

    narrow_transσ_Δμ = narrow_transσ_result.samples.off.μ - narrow_transσ_result.samples.on.μ
    narrow_transσ_title = "Δμ = $(round(narrow_transσ_Δμ,digits=2)), σₒₙ=$(round(narrow_transσ_result.samples.on.σ, digits=2))"
    fig[2,1] = narrow_transσ_ax = Makie.Axis(fig, title=narrow_transσ_title)
    @info """$narrow_transσ_title
    DoS RMSE: $(Optim.minimum(narrow_transσ_result.fits.my))
    Gaussian RMSE: $(Optim.minimum(narrow_transσ_result.fits.meijer))"""

    noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
    label_a = fig[1,1,TopLeft()] = Label(fig, "A", font=noto_sans_bold, halign=:left, textsize=40)
    label_b = fig[1,2,TopLeft()] = Label(fig, "B", font=noto_sans_bold, halign=:left, textsize=40)
    label_c = fig[2,1,TopLeft()] = Label(fig, "C", font=noto_sans_bold, halign=:left, textsize=40)
    label_d = fig[2,2,TopLeft()] = Label(fig, "D", font=noto_sans_bold, halign=:left, textsize=40)
    three_model_pairs = [
        (:my, THREE_MODELS.my),
        (:meijer, THREE_MODELS.meijer)#,
        #(:kim, THREE_MODELS.kim)
    ]
    plot_example_result!(wide_idσ_ax, wide_idσ_result, three_model_pairs)
    plot_example_result!(wide_transσ_ax, wide_transσ_result, three_model_pairs)
    plot_example_result!(narrow_transσ_ax, narrow_transσ_result, three_model_pairs)
    axs = [wide_idσ_ax, wide_transσ_ax, narrow_transσ_ax]
    # axs .|> ax -> ax.xticklabelrotation[] = π/2
    # Makie.Axis(first(drawn_fig)).xticklabelrotation[] = π/2

    axs .|> ax -> (xlims!(ax,(-2,2)); ax.xticks[] = [-1, 0, 1])
    axs .|> ax -> ylims!(ax,(0.0,1.2))
    axs .|> ax -> ax.yticks[] = [0.0, 0.5, 1.0]

    #Colorbar(fig[1,2], only(drawn_fig[1,1].axis.scene.plots))
    display(drawn_fig)

    save(plotsdir("toy_model_fit.svg"), drawn_fig)
end

false

end;
