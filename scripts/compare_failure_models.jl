using NeuralModels: AxisIndices
using DrWatson
@quickactivate "FailureOfInhibition2021"

using Base.Iterators
using Optim, Statistics, AxisIndices

using CairoMakie

using NeuralModels

include(srcdir("axisarray.jl"))

my_model = difference_of_simple_sigmoids
kim_model = product_of_simple_sigmoids
meijer_model = gaussian
models = (my=my_model, kim=kim_model, meijer=meijer_model)
sigmoidal_param_nt = (
    a_on = 1,
    θ_on = 2,
    a_off = 3,
    θ_off = 4
)
gaussian_param_nt = (
    μ=1,
    σ=2,
    A=3
)

model_param_nts = (
    my=sigmoidal_param_nt,
    kim=sigmoidal_param_nt,
    meijer=gaussian_param_nt
)

# source: https://www.wolframcloud.com/env/ecd2d9d4-3c74-44d0-b066-760406ca3f72
# personal mathematica notebook finding f'(μ + sqrt(2log(2))σ) 
# where f is the gaussian PDF, and the argument comes from the FWHM source above
# and a = slope * 4 comes from g'(θ) where g is the sigmoid
function a_from_σ(σ)
    slope = sqrt(log(2) / pi) / 2σ^2
    a = slope * 4
    return a
end
function σ_from_a(a)
    slope = a / 4
    sd = sqrt(sqrt(log(2) / pi) / 2slope)
    return sd
end
function gaussian_get(g, key)
    if key ∈ [:a_on, :a_off]
        a_from_σ(g[gaussian_param_nt.σ])
    elseif key ∈ [:θ_on, :θ_off]
        # source: https://en.wikipedia.org/wiki/Full_width_at_half_maximum
        on_subtract = key == :θ_on ? -1 : 1
        on_subtract * sqrt(2log(2)) * gaussian_get(g, :σ) + gaussian_get(g, :μ)
    else
        g[gaussian_param_nt[key]]
    end
end

function sigmoidal_get(s, key)
    s[sigmoidal_param_nt[key]]
end

function model_get(model_name, m, key)
    if model_name ∈ [:my, :kim]
        sigmoidal_get(m, key)
    elseif model_name == :meijer
        gaussian_get(m, key)
    else
        error("Unknown model $model_name")
    end
end


# NOTE determine whether it makes more sense to link the distributions

randn_mean_sd(T, μ, σ) = randn(T) * σ + μ
fit_error(model_fit, data) = sqrt(mean(@. (model_fit - data)^2))

# Generate off-on-off neurons with normally distributed thresholds, conditional on θ_on < θ_off. Returns (nt, n_discarded) where nt is a namedtuple of fits, and n_discarded says how many random neuron choices were discarded because θ_on >= θ_off
function fit_normally_distributed_thresholds(
        θ_on_mean::T, θ_on_sd,
        θ_off_mean, θ_off_sd, 
        N_neurons=1000, N_xs=1000,
        N_sds_considered=5
    ) where T
    @assert θ_on_mean < θ_off_mean
    lower_bound = min(θ_on_mean - N_sds_considered * θ_on_sd, θ_off_mean - N_sds_considered * θ_off_sd)
    upper_bound = max(θ_on_mean + N_sds_considered * θ_on_sd, θ_off_mean + N_sds_considered * θ_off_sd)
    xs = range(lower_bound, upper_bound, length=N_xs)

    # FIXME no need to compute for all N_xs X N_neurons
    accumulated_N_active = zeros(T, N_xs)
    threshold_on_arr = zeros(T, N_neurons)
    threshold_off_arr = zeros(T, N_neurons)
    accumulated_N_discarded = 0
    for i_neuron in 1:N_neurons
        threshold_on = randn_mean_sd(T, θ_on_mean, θ_on_sd)
        threshold_off = randn_mean_sd(T, θ_off_mean, θ_off_sd)
        while threshold_off <= threshold_on
            accumulated_N_discarded += 1
            threshold_on = randn_mean_sd(T, θ_on_mean, θ_on_sd)
            threshold_off = randn_mean_sd(T, θ_off_mean, θ_off_sd)
        end
        threshold_on_arr[i_neuron] = threshold_on
        threshold_off_arr[i_neuron] = threshold_off
        accumulated_N_active .+= binary_switch_off_on_off.(xs, threshold_on, threshold_off)
    end
    proportion_active = accumulated_N_active ./ N_neurons
    sample_μ_on = mean(threshold_on_arr)
    sample_μ_off = mean(threshold_off_arr)
    sample_σ_on = std(threshold_on_arr)
    sample_σ_off = std(threshold_off_arr)
    #lineplot(xs, proportion_active) |> display

    # Estimate parameters for fitting starting points
    max_proportion_active = maximum(proportion_active)
    guess_θ_on = xs[findfirst(proportion_active .>= max_proportion_active / 2)]
    guess_θ_off = xs[findlast(proportion_active .>= max_proportion_active / 2)]
    guess_plateau_start = xs[findfirst(proportion_active .>= 0.99max_proportion_active)]
    guess_plateau_stop = xs[findlast(proportion_active .>= 0.99max_proportion_active)]
    guess_a_on = 1 / (guess_plateau_start - guess_θ_on)
    guess_a_off = 1 / (guess_θ_off - guess_plateau_stop)
    guess_σ = ((guess_θ_off - guess_θ_on) / 2) / (sqrt(2log(2)))
    guess_μ = guess_θ_on + (guess_θ_off - guess_θ_on) / 2
    guess_A = max_proportion_active

    my_model_loss(p) = fit_error(my_model.(xs, p...), proportion_active)
    kim_model_loss(p) = fit_error(kim_model.(xs, p...), proportion_active)
    meijer_model_loss(p) = fit_error(meijer_model.(xs, p...), proportion_active)

    my_model_guess = [guess_a_on, guess_θ_on, guess_a_off, guess_θ_off]
    kim_model_guess = copy(my_model_guess)
    meijer_model_guess = [guess_μ, guess_σ, guess_A]

    my_model_fit = optimize(my_model_loss, my_model_guess, LBFGS(); autodiff=:forward)
    kim_model_fit = optimize(kim_model_loss, kim_model_guess, LBFGS(); autodiff=:forward)
    meijer_model_fit = optimize(meijer_model_loss, meijer_model_guess, NelderMead())

    return (
        (my=my_model_fit, kim=kim_model_fit, meijer=meijer_model_fit), 
        xs,
        accumulated_N_discarded, 
        (μ_on=sample_μ_on, σ_on=sample_σ_on, μ_off=sample_μ_off, σ_off=sample_σ_off)
    )
end

function fit_range_of_thresholds(
        θ_on_mean, θ_on_sd,
        θ_off_mean, θ_off_sd; 
        kwargs...
    )
    all_axes = (θ_on_mean, θ_on_sd,
                θ_off_mean, θ_off_sd)
    all_axes_names = (:θ_on_mean, :θ_on_sd,
                      :θ_off_mean, :θ_off_sd)
    axes = filter(x -> typeof(x) <: AbstractArray, all_axes)
    axes_names = all_axes_names[(typeof.(all_axes) .<: AbstractArray) |> collect]
    off_on_off_parameters = product(
        θ_on_mean, θ_on_sd,
        θ_off_mean, θ_off_sd
    )
    results = NamedAxisArray{axes_names}(
        map(off_on_off_parameters) do p
            fit_normally_distributed_thresholds(
                p...; kwargs...
            )
        end, 
        axes
            
    )

    return results

end




function vector_plot_models(results::AbstractVector, model_names; title="")
    vals = only(axes_keys(results))
    xlabel = string(only(AxisIndices.dimnames(results)))
    fig = Figure()
    fig[1,1] = minimum_ax = Makie.Axis(fig; title="mse", xlabel=xlabel)
    fig[2,1] = theta_on_ax = Makie.Axis(fig; title="|θ̂_on - μ̂_on|")
    fig[2,2] = theta_off_ax = Makie.Axis(fig; title="|θ̂_off - μ̂_off|")
    fig[3,1] = a_on_ax = Makie.Axis(fig; title="|â_on - a(σ̂_on)|")
    fig[3,2] = a_off_ax = Makie.Axis(fig; title="|â_off - a(σ̂_off)|")
    # fig[3,1] = σ_on_ax = Makie.Axis(fig; title="σ̂_on")
    # fig[3,2] = σ_off_ax = Makie.Axis(fig; title="σ̂_off")
    fig[4,1] = first_example_ax = Makie.Axis(fig; title="first example")
    fig[4,2] = last_example_ax = Makie.Axis(fig; title="last example")
    for name in filter(x -> x != :kim, model_names)
        this_model_results = results .|> x -> x[1][name]
        this_model_samples = results .|> x -> x[4]
        lines!(minimum_ax, vals, Optim.minimum.(this_model_results))
        estimated_params =  Optim.minimizer.(this_model_results)
        θ_on_errs = abs.(
            model_get.(name,estimated_params,:θ_on) .-
            (this_model_samples .|> sample -> sample.μ_on)
        )
        lines!(theta_on_ax, vals, θ_on_errs)
        θ_off_errs = abs.(
            model_get.(name,estimated_params,:θ_off) .-
            (this_model_samples .|> sample -> sample.μ_off)
        )
        lines!(theta_off_ax, vals, θ_off_errs)
        a_on_errs = abs.(
            model_get.(name,estimated_params,:a_on) .-
            (this_model_samples .|> sample -> a_from_σ(sample.σ_on))
        )
        lines!(a_on_ax, vals, a_on_errs)
        a_off_errs = abs.(
            model_get.(name,estimated_params,:a_off) .-
            (this_model_samples .|> sample -> a_from_σ(sample.σ_off))
        )
        lines!(a_off_ax, vals, a_off_errs)
        lines!(first_example_ax, results[begin][2], models[name].(results[begin][2], estimated_params[begin]...))
        lines!(last_example_ax, results[end][2], models[name].(results[end][2], estimated_params[end]...))
    end
    supertitle = fig[0,:] = Label(fig, title)
    display(fig)
end


let θ_on_mean=0., θ_on_sd=1 ./ (2.:10.:52.), 
    θ_off_mean=0.05:0.05:1., θ_off_sd=1 ./ (2.:10.:52.);

results = fit_range_of_thresholds(
    θ_on_mean, θ_on_sd,
    θ_off_mean, θ_off_sd
)

model_names = typeof(results[begin][1]).parameters[1]
high_gain_results = results[θ_off_sd=θ_off_sd[end], θ_on_sd=θ_on_sd[end]]
low_gain_results = results[θ_off_sd=θ_off_sd[begin], θ_on_sd=θ_on_sd[begin]]
wide_hg_off_results = results[θ_off_mean=θ_off_mean[end], θ_off_sd=θ_off_sd[end]]
narrow_hg_off_results = results[θ_off_mean=θ_off_mean[begin], θ_off_sd=θ_off_sd[end]] # varies θ_on_sd


vector_plot_models(high_gain_results, model_names; title="High gain both: off $(θ_off_sd[end]), on:$(θ_on_sd[end])")
vector_plot_models(low_gain_results, model_names; title="Low gain both: off $(θ_off_sd[begin]), on $(θ_on_sd[begin])")
vector_plot_models(wide_hg_off_results, model_names; title="Wide HG Off: θ_off $(θ_off_mean[end]), σ_off $(θ_off_sd[end])")
vector_plot_models(narrow_hg_off_results, model_names; title="Narrow HG Off: θ_off $(θ_off_mean[begin]), σ_off $(θ_off_sd[end])")

end
