using Optim, Statistics, AxisIndices

using CairoMakie,Makie

using NeuralModels

using Memoize

include(srcdir("axisarray.jl"))
include(srcdir("estimate_params.jl"))
include(srcdir("theme.jl"))

if !@isdefined(force_results_calc) || force_results_calc
    μ_on=0.; σ_on=0.04:0.04:0.36; 
    μ_off=0.1:0.1:1.; σ_off=0.04:0.04:0.36    
    results = fit_range_of_thresholds(;
        μ_on=μ_on, σ_on=σ_on,
        μ_off=μ_off, σ_off=σ_off
    )
    force_results_calc = false
end

force_results_calc = let results = results;
model_string_names = Dict(
    :my => "DoS",
    :meijer => "Gauss",
    :kim => "Kim"
)
model_rmses = Dict(
    key => [Optim.minimum(getproperty(result.fits, key)) for result in results] for key in keys(model_string_names)
)
model_rmses_means = Dict(
    key => mean(model_rmses[key]) for key in keys(model_rmses)
)
model_rmses_sds = Dict(
    key => std(model_rmses[key]) for key in keys(model_rmses)
)

for key in keys(model_rmses)
    @info "$(model_string_names[key]) RMSE: $(model_rmses_means[key]) ± $(model_rmses_sds[key])"
end

results_df = results_nameddimsarray_to_df(results; model_names=[:my, :meijer, :kim])

sample_μ_vs_μ_plt = data(results_df) * mapping(
    :μ_off,
    (:samples) => ((s) -> s.off.μ - s.on.μ) => "Δμ",
) * expectation() + data(results_df) * mapping(:μ_off, :μ_off => "Δμ") * expectation() * visual(Lines)
display(draw(sample_μ_vs_μ_plt))

sample_μ_vs_diffs_plt = data(results_df) * mapping(
    :μ_off,
    (:σ_on, :σ_off) => ((son, soff) -> abs(son - soff)) => "Δσ",
    (:samples) => ((s) -> s.off.μ - s.on.μ) => "Δμ",
) * expectation()
display(draw(sample_μ_vs_diffs_plt))

sample_μ_err_vs_μ_plt = data(results_df) * mapping(
    :μ_off,
    (:samples, :μ_off) => ((s, μ) -> s.off.μ - s.on.μ - μ) => "Δμ error",
) * expectation()
display(draw(sample_μ_err_vs_μ_plt))

sample_μ_err_vs_diffs_plt = data(results_df) * mapping(
    :μ_off,
    (:σ_on, :σ_off) => ((son, soff) -> abs(son - soff)) => "Δσ",
    (:samples, :μ_off) => ((s, μ) -> s.off.μ - s.on.μ - μ) => "Δμ error",
) * expectation()
display(draw(sample_μ_err_vs_diffs_plt))

n_discarded_vs_diffs_plt = data(results_df) * mapping(
    :μ_off,
    (:σ_on, :σ_off) => ((son, soff) -> abs(son - soff)) => "Δσ",
    (:N_discarded) => "discarded neurons",
) * expectation()
display(draw(n_discarded_vs_diffs_plt))

false
end;