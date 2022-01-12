
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