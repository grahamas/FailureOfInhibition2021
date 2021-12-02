
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