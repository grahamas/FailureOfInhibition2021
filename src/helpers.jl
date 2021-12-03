
using NamedDims

function avg_across_dims(arr, dims)
    avgd = mean_skip_missing(arr; dims=dims)
    squeezed = dropdims(avgd, dims=dims)
    return squeezed
end

function _collapse_to_axes(A::NamedDimsArray{name_syms}, x_sym::Symbol, y_sym::Symbol) where name_syms
    collapsed_syms = Tuple(setdiff(name_syms, (y_sym, x_sym)))
    data = if findfirst(name_syms .== y_sym) < findfirst(name_syms .== x_sym)
        Simulation73Plotting.avg_across_dims(A, collapsed_syms)'
    else
        Simulation73Plotting.avg_across_dims(A, collapsed_syms)
    end
    return data
end

function _collapse_to_axes(A::NamedDimsArray{name_syms}, x_sym::Symbol) where name_syms
    collapsed_syms = Tuple(setdiff(name_syms, [x_sym]))
    return Simulation73Plotting.avg_across_dims(A, collapsed_syms)
end

function enumerate_nda_dims(nda::NamedDimsArray{NAMES}, dims::NamedTuple{NAMES}) where NAMES
    zip(NamedTuple{NAMES}.(product(dims...)), nda)
end


function subset_axes(axes_nt::NamedTuple{Names}, subset_range) where Names
    NamedTuple{Names}(
        [axis[axis .∈ subset_range] for axis in axes_nt]
    )
end

function subset_nda_dims(nda::NamedDimsArray{Names}, nda_dims::NamedTuple{Names}, subset_range) where Names
    logical_subsets = NamedTuple{Names}(
        [axis .∈ subset_range for axis in nda_dims]
    )
    subset_nda = getindex(nda; logical_subsets...)
    subset_nda_dims = NamedTuple{Names}(
        [nda_dims[name][logical_subsets[name]] for name in Names]
    )
    return (subset_nda, subset_nda_dims)
end

function abbrev_count_label(x)
    if x >= 1000
        try
            "$(Int(x / 1000))K"
        catch
            "$(x / 1000)K"
        end
    else
        "$(Int(x))"
    end
end

function max_percentile_conditions(values; left_tail_pct=0.025, right_tail_pct=0.025)
    # Note: not symmetric tails, since values can be negative
    left_tail = maximum(values) * left_tail_pct
    right_tail = maximum(values) * (1 - right_tail_pct)
    middle_pct = 1 - (left_tail_pct + right_tail_pct)
    conditions = map(values) do value
        if value <= left_tail
            "min"
        elseif left_tail < value < right_tail
            "mid"
        elseif right_tail <= value
            "max"
        else
            throw(DomainError(value))
        end
    end
    return (conditions, left_tail, right_tail)
end

format_percent(x) = "$(round(Int, 100x))%"
format_tail(x) = "$(round(x, sigdigits=2))"