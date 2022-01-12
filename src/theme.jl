using ColorSchemes, ColorTypes

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
    fontsize=56,
    strokewidth= 5.,
    Axis = (
        backgroundcolor = RGBA(1.,1.,1.,0.),
        leftspinevisible = true,
        rightspinevisible = false,
        bottomspinevisible = true,
        topspinevisible = false,
        xgridcolor = RGBA(1.,1.,1.,0.),
        ygridcolor = RGBA(1.,1.,1.,0.),
        strokewidth= 5.,
        ytickformat = xs -> abbrev_count_label.(xs)
    )
)

nullcline_theme = Theme(
    fontsize=48,
    Axis = (
        backgroundcolor = RGBA(1.,1.,1.,0.),
        leftspinevisible = true,
        rightspinevisible = false,
        bottomspinevisible = true,
        topspinevisible = false,
        strokewidth=2.,
        xgridcolor = RGBA(1.,1.,1.,0.),
        ygridcolor = RGBA(1.,1.,1.,0.)
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

violin_theme = Theme(
    fontsize=24,
    Axis = (
        backgroundcolor = RGBA(1.,1.,1.,0.),
        leftspinevisible = true,
        rightspinevisible = false,
        bottomspinevisible = true,
        topspinevisible = false,
        strokewidth=2.,
        xgridcolor = RGBA(1.,1.,1.,0.),
        ygridcolor = RGBA(1.,1.,1.,0.)
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
    ),
    Violin = (
        bandwidth=0.003,
        npoints=256*16
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