using Base: sign_mask
using Parameters
using AlgebraOfGraphics, DataFrames

abstract type NonlP{T} end
@with_kw struct TwoSigmoidsP{T} <: NonlP{T}
    a_on::T
    θ_on::T
    a_off::T
    θ_off::T
end
Base.collect(p::TwoSigmoidsP) = [p.a_on, p.θ_on, p.a_off, p.θ_off]
TwoSigmoidsP(vec::Vector{T}) where T = TwoSigmoidsP{T}(vec...)
@with_kw struct GaussianP{T} <: NonlP{T}
    μ::T
    σ::T
    A::Union{T,Missing}
end
Base.collect(p::GaussianP) = [p.μ, p.σ, p.A]
GaussianP(vec::Vector{T}) where T = GaussianP{T}(vec...)

# source: https://www.wolframcloud.com/env/ecd2d9d4-3c74-44d0-b066-760406ca3f72
# personal mathematica notebook finding f'(μ + sqrt(2log(2))σ) 
# where f is the gaussian PDF, and the argument comes from the FWHM source above
# and a = slope * 4 comes from g'(θ) where g is the sigmoid
# THIS ONLY APPLIES TO FITTING GAUSSIAN TO WHOLE UP-DOWN
function a_from_gaussian_fit_σ(σ)
    slope = sqrt(log(2) / pi) / 2σ^2
    a = slope * 4
    return a
end
function gaussian_fit_σ_from_a(a)
    slope = a / 4
    sd = sqrt(sqrt(log(2) / pi) / 2slope)
    return sd
end
function Base.getproperty(g::GaussianP, key::Symbol)
    if key ∈ [:a_on, :a_off]
        a_from_gaussian_fit_σ(getfield(g, :σ))
    elseif key ∈ [:θ_on, :θ_off, :μ_on, :μ_off]
        # source: https://en.wikipedia.org/wiki/Full_width_at_half_maximum
        on_subtract = key ∈ [:θ_on, :μ_on] ? -1 : 1
        on_subtract * sqrt(2log(2)) * getfield(g, :σ) + getfield(g, :μ)
    elseif key ∈ [:σ_on, :σ_off]
        getfield(g, :σ)
    else
        getfield(g, key)
    end
end
function Base.getproperty(s::TwoSigmoidsP, key::Symbol)
    if key == :σ_on
        sample_σ_from_a(getfield(s, :a_on))
    elseif key == :σ_off
        sample_σ_from_a(getfield(s, :a_off))
    elseif key == :μ_on
        getfield(s, :θ_on)
    elseif key == :μ_off
        getfield(s, :θ_off)
    else
        getfield(s, key)
    end
end

# max slope of sigmoid == peak of gaussian
# sigmoid `a` param == slope * 4
function a_from_sample_σ(σ)
    max_slope = 1 / (σ * sqrt(2π))
    a = max_slope * 4
    return a
end

function sample_σ_from_a(a)
    max_slope = a / 4
    σ = 1 / (max_slope * sqrt(2π)) 
    return σ
end


# NOTE determine whether it makes more sense to link the distributions
randn_mean_sd(T, μ, σ) = (randn(T) * σ) + μ
fit_error(model_fit, data) = sqrt(mean(@. (model_fit - data)^2))

# Generate off-on-off neurons with normally distributed thresholds, conditional on θ_on < θ_off. Returns (nt, n_discarded) where nt is a namedtuple of fits, and n_discarded says how many random neuron choices were discarded because θ_on >= θ_off
function fit_normally_distributed_thresholds(
        μ_on::T, σ_on,
        μ_off, σ_off, 
        N_neurons=1000, N_xs=1000,
        N_sds_considered=5
    ) where T
    @assert μ_on < μ_off
    lower_bound = min(μ_on - N_sds_considered * σ_on, μ_off - N_sds_considered * σ_off)
    upper_bound = max(μ_on + N_sds_considered * σ_on, μ_off + N_sds_considered * σ_off)
    xs = range(lower_bound, upper_bound, length=N_xs)

    # FIXME no need to compute for all N_xs X N_neurons
    accumulated_N_active = zeros(T, N_xs)
    threshold_on_arr = zeros(T, N_neurons)
    threshold_off_arr = zeros(T, N_neurons)
    accumulated_N_discarded = 0
    for i_neuron in 1:N_neurons
        threshold_on = randn_mean_sd(T, μ_on, σ_on)
        threshold_off = randn_mean_sd(T, μ_off, σ_off)
        while threshold_off <= threshold_on
            accumulated_N_discarded += 1
            threshold_on = randn_mean_sd(T, μ_on, σ_on)
            threshold_off = randn_mean_sd(T, μ_off, σ_off)
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
    guess_a_on = 1/ (guess_plateau_start - guess_θ_on)
    guess_a_off = 1/(guess_θ_off - guess_plateau_stop)
    guess_σ = ((guess_θ_off - guess_θ_on) / 2) / (sqrt(2log(2)))
    guess_μ = guess_θ_on + (guess_θ_off - guess_θ_on) / 2
    guess_A = max_proportion_active

    my_model_loss(p) = fit_error(
        my_model.(xs,Ref(TwoSigmoidsP(p))),
        proportion_active
    )
    kim_model_loss(p) = fit_error(
        kim_model.(xs,Ref(TwoSigmoidsP(p))),
        proportion_active
    )
    meijer_model_loss(p) = fit_error(
        meijer_model.(xs,Ref(GaussianP(p))),
        proportion_active
    )

    my_model_guess = TwoSigmoidsP(;
        a_on=guess_a_on, θ_on=guess_θ_on, 
        a_off=guess_a_off, θ_off=guess_θ_off
    ) |> collect
    kim_model_guess = TwoSigmoidsP(;
        a_on=guess_a_on, θ_on=guess_θ_on, 
        a_off=guess_a_off, θ_off=guess_θ_off
    ) |> collect
    meijer_model_guess = GaussianP(;
        μ=guess_μ, σ=guess_σ, A=guess_A
    ) |> collect

    my_model_fit = optimize(my_model_loss, my_model_guess, LBFGS(); autodiff=:forward)
    kim_model_fit = optimize(kim_model_loss, kim_model_guess, LBFGS(); autodiff=:forward)
    meijer_model_fit = optimize(meijer_model_loss, meijer_model_guess, NelderMead())

    return (
        fits=(
            my=my_model_fit,
            kim=kim_model_fit,
            meijer=meijer_model_fit
        ),
        minimizing_p=(
            my=TwoSigmoidsP(Optim.minimizer(my_model_fit)), 
            kim=TwoSigmoidsP(Optim.minimizer(kim_model_fit)), 
            meijer=GaussianP(Optim.minimizer(meijer_model_fit))
        ),
        xs=xs,
        N_discarded=accumulated_N_discarded, 
        samples=(
            on=GaussianP(μ=sample_μ_on, σ=sample_σ_on, A=missing),
            off=GaussianP(μ=sample_μ_off, σ=σ_off=sample_σ_off, A=missing)
        ),
        truth=proportion_active
    )
end

@memoize function fit_range_of_thresholds(;
        μ_on, σ_on,
        μ_off, σ_off, 
        kwargs...
    )
    all_axes = (μ_on, σ_on,
                μ_off, σ_off)
    all_axes_names = (:μ_on, :σ_on,
                      :μ_off, :σ_off)
    axes = filter(x -> typeof(x) <: AbstractArray, all_axes)
    axes_names = all_axes_names[(typeof.(all_axes) .<: AbstractArray) |> collect]
    off_on_off_parameters = product(
        μ_on, σ_on,
        μ_off, σ_off
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
    # fig[3,1] = a_on_ax = Makie.Axis(fig; title="|â_on - a(σ̂_on)|")
    # fig[3,2] = a_off_ax = Makie.Axis(fig; title="|â_off - a(σ̂_off)|")
    fig[3,1] = σ_on_ax = Makie.Axis(fig; title="|σ̂_on - σ(â_on)|")
    fig[3,2] = σ_off_ax = Makie.Axis(fig; title="|σ̂_off - σ(â_off)|")
    fig[4,1] = first_example_ax = Makie.Axis(fig; title="first example: $(vals[begin])")
    fig[4,2] = last_example_ax = Makie.Axis(fig; title="last example: $(vals[end])")
    plotted_names = filter(x -> x != :kim, model_names) |> collect
    last_example_lines = map(plotted_names) do name
        minima = results .|> x -> Optim.minimum(x.fits[name])
        this_model_samples = results .|> x -> x.samples
        lines!(minimum_ax, vals, minima)
        estimated_params = results .|> result -> result.minimizing_p[name]
        θ_on_errs = abs.(
            getproperty.(estimated_params, :θ_on) .-
            (this_model_samples .|> sample -> sample.on.μ)
        )
        lines!(theta_on_ax, vals, θ_on_errs)
        θ_off_errs = abs.(
            getproperty.(estimated_params, :θ_off) .-
            (this_model_samples .|> sample -> sample.off.μ)
        )
        lines!(theta_off_ax, vals, θ_off_errs)
        # a_on_errs = abs.(
        #     getproperty.(estimated_params, :a_on) .-
        #     (this_model_samples .|> sample -> abs(a_from_sample_σ(sample.on.σ)))
        # )
        # lines!(a_on_ax, vals, a_on_errs)
        # a_off_errs = abs.(
        #     getproperty.(estimated_params, :a_off) .-
        #     (this_model_samples .|> sample -> abs(a_from_sample_σ(sample.off.σ)))
        # )
        # lines!(a_off_ax, vals, a_off_errs)
        σ_on_errs = abs.(
            getproperty.(estimated_params, :σ_on) .-
            (this_model_samples .|> sample -> abs(sample.on.σ))
        )
        lines!(σ_on_ax, vals, σ_on_errs)
        σ_off_errs = abs.(
            getproperty.(estimated_params, :σ_off) .-
            (this_model_samples .|> sample -> abs(sample.off.σ))
        )
        lines!(σ_off_ax, vals, σ_off_errs)
        lines!(first_example_ax, results[begin].xs, models[name].(results[begin].xs, Ref(estimated_params[begin])))
        lines!(last_example_ax, results[end].xs, models[name].(results[end].xs, Ref(estimated_params[end])))
    end
    lines!(first_example_ax, results[begin].xs, results[begin].truth)
    true_last_example_line = lines!(last_example_ax, results[end].xs, results[end].truth)
    Legend(fig[1,2], [last_example_lines..., true_last_example_line], [string.(plotted_names)..., "truth"], tellwidth=false, tellheight=false)
    supertitle = fig[0,:] = Label(fig, title)
    display(fig)
end

nt_get(nt, keys::Union{Tuple,AbstractVector}) = NamedTuple{keys}([nt[key] for key in keys])

function surplus_error_df(arr::NamedAxisArray, 
        baseline_est_fn, baseline_truth_fn,
        surplus_est_fn, surplus_truth_fn,
        surplus_error_name=:surplus_error)
    baseline_est_arr = baseline_est_fn.(arr)
    @show maximum(baseline_est_arr)
    baseline_truth_arr = baseline_truth_fn.(arr)
    surplus_est_arr = surplus_est_fn.(arr)
    @show maximum(surplus_est_arr)
    surplus_truth_arr = surplus_truth_fn.(arr)
    surplus_error_arr = @. abs(surplus_est_arr - surplus_truth_arr) - abs(baseline_est_arr - baseline_truth_arr)
    return axisarray_to_df(surplus_error_arr; value_col_name=surplus_error_name)
end

function axisarray_to_df(arr::NamedAxisArray{NAMES,T}; value_col_name::Symbol=:value) where {NAMES,T}
    first_coord, first_val = first(enumerate_naa(arr))
    N_coords = length(first_coord)
    @show T
    @show typeof(first_coord)
    @show typeof(arr)
    @assert all(eltype(first_coord).parameters .== T)
    coord_col_names = keys(first_coord)
    df_arr = Matrix{T}(undef, length(arr), N_coords + 1)
    for (i, (coord, val)) in enumerate(enumerate_naa(arr))
        df_arr[i, 1:N_coords] .= collect(coord)
        df_arr[i, N_coords+1] = val
    end
    df = DataFrame(; Dict(coord_col_names[i] => df_arr[:,i] for i ∈ 1:length(coord_col_names))..., Dict(value_col_name => df_arr[:, N_coords+1])...)
    return df
end

@memoize function results_axisarray_to_df(arr::NamedAxisArray; title="", model_names=nothing)
    first_coord, first_val = first(enumerate_naa(arr))
    @assert first_val isa NamedTuple
    val_names = keys(first_val)
    @assert val_names == (:fits, :minimizing_p, :xs, :N_discarded, :samples, :truth)
    model_variables = (:fits, :minimizing_p)
    model_invariants = filter(name -> name ∉ model_variables, val_names)
    val_types = map(val_names) do name
        if name ∈ model_variables
            Union{typeof.(collect(first_val[name]))...}
        else
            typeof(first_val[name])
        end
    end
    model_names = model_names === nothing ? keys(first_val.fits) : model_names
    df_val_fields = map(zip(val_names, val_types)) do (name, type)
        name => type[]
    end
    df_coord_fields = map(zip(keys(first_coord), typeof.(values(first_coord)))) do (name, type)
        name => type[]
    end
    df = DataFrame(; df_coord_fields..., df_val_fields..., model=Symbol[])
    for (coord, val) in enumerate_naa(arr)
        same_across_models = nt_get(val, model_invariants)
        for model_name in model_names
            different_across_models = NamedTuple{model_variables}(
                [val[var][model_name] for var in model_variables]
            )
            push!(df, reduce(merge, [coord, same_across_models, different_across_models, Dict(:model => model_name)]))
        end
    end
    return df
end

function fit_vs_differences(results::AbstractArray; 
        model_names, kwargs...
    )
    df = results_axisarray_to_df(results; model_names=model_names)
    fit_vs_differences(df; kwargs...)
end
function fit_vs_differences(df::DataFrame; 
    title="", μ_on
)
    df.μ_difference = abs.(df.μ_off .- μ_on)
    df.σ_difference = abs.(df.σ_off .- df.σ_on)
    df.mse = Optim.minimum.(df.fits)
    mse_plt = data(df) * mapping(:μ_difference, :σ_difference, :mse, layout=:model) * expectation()
    df.μ_on_error = abs.(
            (df.minimizing_p .|> p -> p.μ_on) .-
            (df.samples .|> s -> s.on.μ)
    )
    μ_on_plt = data(df) * mapping(:μ_difference, :σ_difference, :μ_on_error, layout=:model) * expectation()
    df.μ_off_error = abs.(
            (df.minimizing_p .|> p -> p.μ_off) .-
            (df.samples .|> s -> s.off.μ)
    )
    μ_off_plt = data(df) * mapping(:μ_difference, :σ_difference, :μ_off_error, layout=:model) * expectation()
    df.σ_on_error = abs.(
            (df.minimizing_p .|> p -> p.σ_on) .-
            (df.samples .|> s -> s.on.σ)
    )
    σ_on_plt = data(df) * mapping(:μ_difference, :σ_difference, :σ_on_error, layout=:model) * expectation()
    df.σ_off_error = abs.(
            (df.minimizing_p .|> p -> p.σ_off) .-
            (df.samples .|> s -> s.off.σ)
    )
    σ_off_plt = data(df) * mapping(:μ_difference, :σ_difference, :σ_off_error, layout=:model) * expectation()
    return [mse_plt, μ_on_plt, μ_off_plt, σ_on_plt, σ_off_plt]
end

function err_visual(arr)
    visual(colorrange=(min(0., minimum(arr)), maximum(arr)))
end

function surplus_error_vs_differences(results::AbstractArray, μ_on)
    axis_labels = (
        :μ_difference => "|μ(off) - μ(on)|", 
        :σ_difference => "|σ(off) - σ(on)|" 
    )

    # Surplus MSE
    surplus_mse_df = surplus_error_df(results, 
        res -> Optim.minimum(res.fits.my), res -> 0.,
        res -> Optim.minimum(res.fits.meijer), res -> 0.,
        :surplus_mse
    )
    surplus_mse_df.μ_difference = abs.(surplus_mse_df.μ_off .- μ_on)
    surplus_mse_df.σ_difference = abs.(surplus_mse_df.σ_off .- surplus_mse_df.σ_on)
    surplus_mse_plot = data(surplus_mse_df) * mapping(axis_labels..., :surplus_mse) * expectation() * err_visual(surplus_mse_df.surplus_mse)
    # Surplus σ_off
    surplus_σ_off_error_df = surplus_error_df(results, 
        res -> res.minimizing_p.my.σ_off, res -> res.samples.off.σ,
        res -> res.minimizing_p.meijer.σ_off, res -> res.samples.off.σ,
        :surplus_σ_off_error
    )
    surplus_σ_off_error_df.μ_difference = abs.(surplus_σ_off_error_df.μ_off .- μ_on)
    surplus_σ_off_error_df.σ_difference = abs.(surplus_σ_off_error_df.σ_off .- surplus_σ_off_error_df.σ_on)
    surplus_σ_off_plot = data(surplus_σ_off_error_df) * mapping(axis_labels..., :surplus_σ_off_error) * expectation() * err_visual(surplus_σ_off_error_df.surplus_σ_off_error)

    # Surplus σ_on
    surplus_σ_on_error_df = surplus_error_df(results, 
        res -> res.minimizing_p.my.σ_on, res -> res.samples.on.σ,
        res -> res.minimizing_p.meijer.σ_on, res -> res.samples.on.σ,
        :surplus_σ_on_error
    )
    surplus_σ_on_error_df.μ_difference = abs.(surplus_σ_on_error_df.μ_off .- μ_on)
    surplus_σ_on_error_df.σ_difference = abs.(surplus_σ_on_error_df.σ_off .- surplus_σ_on_error_df.σ_on)
    surplus_σ_on_plot = data(surplus_σ_on_error_df) * mapping(axis_labels..., :surplus_σ_on_error) * expectation() * err_visual(surplus_σ_on_error_df.surplus_σ_on_error)

    # Surplus μ_off
    surplus_μ_off_error_df = surplus_error_df(results, 
        res -> res.minimizing_p.my.μ_off, res -> res.samples.off.μ,
        res -> res.minimizing_p.meijer.μ_off, res -> res.samples.off.μ,
        :surplus_μ_off_error
    )
    surplus_μ_off_error_df.μ_difference = abs.(surplus_μ_off_error_df.μ_off .- μ_on)
    surplus_μ_off_error_df.σ_difference = abs.(surplus_μ_off_error_df.σ_off .- surplus_μ_off_error_df.σ_on)
    surplus_μ_off_plot = data(surplus_μ_off_error_df) * mapping(axis_labels..., :surplus_μ_off_error) * expectation() * err_visual(surplus_μ_off_error_df.surplus_μ_off_error)

    # Surplus μ_on
    surplus_μ_on_error_df = surplus_error_df(results, 
        res -> res.minimizing_p.my.μ_on, res -> res.samples.on.μ,
        res -> res.minimizing_p.meijer.μ_on, res -> res.samples.on.μ,
        :surplus_μ_on_error
    )
    surplus_μ_on_error_df.μ_difference = abs.(surplus_μ_on_error_df.μ_off .- μ_on)
    surplus_μ_on_error_df.σ_difference = abs.(surplus_μ_on_error_df.σ_off .- surplus_μ_on_error_df.σ_on)
    surplus_μ_on_plot = data(surplus_μ_on_error_df) * mapping(axis_labels..., :surplus_μ_on_error) * expectation() * err_visual(surplus_μ_on_error_df.surplus_μ_on_error)

    return (surplus_mse=surplus_mse_plot, 
        surplus_μ_on_error=surplus_μ_on_plot,
        surplus_σ_on_error=surplus_σ_on_plot,
        surplus_μ_off_error=surplus_μ_off_plot,
        surplus_σ_off_error=surplus_σ_off_plot)
end
