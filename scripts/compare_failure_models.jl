using DrWatson
@quickactivate "FailureOfInhibition2021"

using UnicodePlots

using NeuralModels

my_model = difference_of_simple_sigmoids
kim_model = product_of_simple_sigmoids
meijer_model = gaussian

# NOTE determine whether it makes more sense to link the distributions

randn_mean_sd(T, μ, σ) = randn(T) * σ + μ
fit_error(model_fit, data) = sqrt(sum(@. (model_fit - data)^2))

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
    accumulated_N_discarded = 0
    for i_neuron in 1:N_neurons
        threshold_on = randn_mean_sd(T, θ_on_mean, θ_on_sd)
        threshold_off = randn_mean_sd(T, θ_off_mean, θ_off_sd)
        while threshold_off <= threshold_on
            accumulated_N_discarded += 1
            threshold_on = randn_mean_sd(T, θ_on_mean, θ_on_sd)
            threshold_off = randn_mean_sd(T, θ_off_mean, θ_off_sd)
        end
        accumulated_N_active .+= binary_switch_off_on_off.(xs, threshold_on, threshold_off)
    end
    proportion_active = accumulated_N_active ./ N_neurons
    lineplot(xs, proportion_active) |> display

    # Estimate parameters for fitting starting points
    max_proportion_active = maximum(proportion_active)
    guess_θ_on = xs[findfirst(proportion_active .>= max_proportion_active / 2)]
    guess_θ_off = xs[findlast(proportion_active .>= max_proportion_active / 2)]
    guess_plateau_start = xs[findfirst(proportion_active .>= 0.99max_proportion_active)]
    guess_plateau_stop = xs[findlast(proportion_active .>= 0.99max_proportion_active)]
    guess_a_on = 1 / (guess_plateau_start - guess_θ_on)
    guess_a_off = 1 / (guess_θ_off - guess_plateau_stop)
    guess_σ = (guess_θ_off - guess_θ_on) / 2
    guess_μ = guess_θ_on + guess_σ

    my_model_loss(p) = fit_error(my_model.(xs, p...), proportion_active)
    kim_model_loss(p) = fit_error(kim_model.(xs, p...), proportion_active)
    meijer_model_loss(p) = fit_error(meijer_model.(xs, p...), proportion_active)

    my_model_guess = [guess_a_on, guess_θ_on, guess_a_off, guess_θ_off]
    kim_model_guess = copy(my_model_guess)
    meijer_model_guess = [guess_μ, guess_σ]

    my_model_fit = optimize(my_model_loss, my_model_guess, LBFGS(); autodiff=:forward)
    kim_model_fit = optimize(kim_model_loss, kim_model_guess, LBFGS(); autodiff=:forward)
    meijer_model_fit = optimize(meijer_model_loss, meijer_model_guess, NelderMead())

    return (my_model_fit, kim_model_fit, meijer_model_fit)
end