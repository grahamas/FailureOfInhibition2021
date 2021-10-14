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
    fontsize=48,
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
        linewidth=7.0,
    ),
    Arrows = (
        arrowsize=20, lengthscale=0.017,
        linewidth=2,
        arrowcolor=:black, linecolor=:black,
        colormap=ColorSchemes.Greys_5,
        normalize=false
    ),
    Scatter = (
        markersize=50,
        strokewidth=2
    ),
    Legend = (
        labelsize = 20,
        textsize = 22,
        linewidth=5.
    )
)  


model_comparison_theme = Theme(
    fontsize=48,
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
        linewidth=9.0,
    ),
    Band = (
        linewidth=20.0,
    ),
    # Arrows = (
    #     arrowsize=20, lengthscale=0.017,
    #     linewidth=2,
    #     arrowcolor=:black, linecolor=:black,
    #     colormap=ColorSchemes.Greys_5,
    #     normalize=false
    # ),
    Scatter = (
        markersize=10.,
        strokewidth=0.
    ),
    Legend = (
        labelsize = 40,
        fontsize = 40,
        linewidth=20.
    ),
    Label = (
        textsize = 56,
    )
)  