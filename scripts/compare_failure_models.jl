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

function my_model(xs,p::TwoSigmoidsP)
    difference_of_simple_sigmoids(xs,p.a_on,p.θ_on,p.a_off,p.θ_off)
end
function kim_model(xs,p::TwoSigmoidsP)
    product_of_simple_sigmoids(xs,p.a_on,p.θ_on,p.a_off,p.θ_off)
end
function meijer_model(xs,p::GaussianP)
    gaussian(xs,p.σ,p.μ,p.A)
end
models = (
    my=my_model,
    kim=kim_model,
    meijer=meijer_model
)

results_df = let μ_on=0., σ_on=sort(1 ./ (1.:3.:25.)), 
    μ_off=0.05:0.05:1., σ_off=sort(1 ./ (1.:3.:25.));

results = fit_range_of_thresholds(;
    μ_on=μ_on, σ_on=σ_on,
    μ_off=μ_off, σ_off=σ_off
)

# model_names = typeof(results[begin][1]).parameters[1]
# high_gain_results = results[σ_off=σ_off[begin], σ_on=σ_on[begin]]
# low_gain_results = results[σ_off=σ_off[end], σ_on=σ_on[end]]
# wide_hg_off_results = results[μ_off=μ_off[end], σ_off=σ_off[begin]]
# narrow_hg_off_results = results[μ_off=μ_off[begin], σ_off=σ_off[begin]] # varies σ_on


# vector_plot_models(high_gain_results, model_names; title="High gain both: off $(σ_off[begin]), on:$(σ_on[begin])")
# vector_plot_models(low_gain_results, model_names; title="Low gain both: off $(σ_off[end]), on $(σ_on[end])")
# vector_plot_models(wide_hg_off_results, model_names; title="Wide HG Off: θ_off $(μ_off[end]), σ_off $(σ_off[begin])")
# vector_plot_models(narrow_hg_off_results, model_names; title="Narrow HG Off: θ_off $(μ_off[begin]), σ_off $(σ_off[begin])")

# df = results_axisarray_to_df(results; model_names=[:my, :meijer])
# plts = fit_vs_differences(df; μ_on=μ_on)
# for i in 1:length(plts)
#     fig = Figure()
#     draw!(fig[1,1], plts[i])
#     display(fig)    
# end

# df
surplus_plots_nt = surplus_error_vs_differences(results, μ_on)

with_theme(nullcline_theme) do
for name in keys(surplus_plots_nt)
    fig = Figure(title=string(name))
    drawn_fig = draw!(fig[1,1], surplus_plots_nt[name]
    #; 
    # fontsize=30, (
    #     backgroundcolor = :white,
    #     leftspinevisible = true,
    #     rightspinevisible = false,
    #     bottomspinevisible = true,
    #     topspinevisible = false,
    #     xgridcolor = :white,
    #     ygridcolor = :white,
    #     ytickformat = xs -> abbrev_count_label.(xs)
    # )...
    )
    Colorbar(fig[1,2], only(drawn_fig[1,1].axis.scene.plots))
    Makie.Axis(first(drawn_fig)).title[] = string(name)
    display(fig)
end
end

surplus_plots_nt

end;
