using ColorSchemes

abbrev_count_label = x -> begin
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

bar_theme = Theme(
    fontsize=30,
    Axis = (
        backgroundcolor = :white,
        leftspinevisible = true,
        rightspinevisible = false,
        bottomspinevisible = true,
        topspinevisible = false,
        xgridcolor = :white,
        ygridcolor = :white,
        ytickformat = xs -> abbrev_count_label.(xs)
    )
)
nullcline_theme = Theme(
    fontsize=30,
    Axis = (
        backgroundcolor = :white,
        leftspinevisible = true,
        rightspinevisible = false,
        bottomspinevisible = true,
        topspinevisible = false,
        xgridcolor = :white,
        ygridcolor = :white
    ),
    Lines = (
        linewidth=4.0,
    ),
    Arrows = (
        arrowsize=10, lengthscale=0.017,
        linewidth=2,
        arrowcolor=:black, linecolor=:black,
        colormap=ColorSchemes.Greys_5,
        normalize=false
    ),
    Scatter = (
        markersize=27,
        strokewidth=1
    )
)  