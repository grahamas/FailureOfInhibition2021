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

# surplus_plots_nt = surplus_error_vs_differences(results, μ_on,
#     [:μ_off => (off -> off - μ_on) => "μ(off) - μ(on)"]
# )

# plot_title_map = (
#     surplus_mse="MSE(Gaussian) - MSE(DoS)", 
#     surplus_μ_on_error="ε(μ_on)(Gaussian) - ε(μ_on)(DoS)",
#     surplus_σ_on_error="ε(σ_on)(Gaussian) - ε(σ_on)(DoS)",
#     surplus_μ_off_error="ε(μ_off)(Gaussian) - ε(μ_off)(DoS)",
#     surplus_σ_off_error="ε(σ_off)(Gaussian) - ε(σ_off)(DoS)"
# )

# with_theme(nullcline_theme) do
# for name in keys(surplus_plots_nt)
#     fig = Figure(title=string(name))
#     drawn_fig = draw!(fig[1,1], surplus_plots_nt[name])
#     #Colorbar(fig[1,2], only(drawn_fig[1,1].axis.scene.plots))
#     Makie.Axis(first(drawn_fig)).title[] = string(plot_title_map[name])
#     Makie.Axis(first(drawn_fig)).xticklabelrotation[] = π/2
#     display(fig)
# end
# end

surplus_plots_nt

end;
