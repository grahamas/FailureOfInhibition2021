using Contour, GeometryTypes

function lifted_wcm_param(;
    Aee=1., Aie=1., Aei=1.5, Aii=0.25,
    θef=0.125, θif=0.4, θib=7., τ=0.4, β=50.,
    decaye=1., decayi=1.)
    @lift WCMMonParams(;
            Aee=$(Node(Aee)), Aie=$(Node(Aie)), Aei=$(Node(Aei)), Aii=$(Node(Aii)),
            θef=$(Node(θef)), θif=$(Node(θif)), θib=$(Node(θib)), τ=$(Node(τ)), β=$(Node(β)),
            decaye=$(Node(decaye)), decayi=$(Node(decayi))
        )
end

function plot_nullclines!(fig::Figure, p::AbstractWCMNullclineParams, axis_length::Integer=100; kwargs...)
    plot_nullclines!(fig, p, [range(-eps(Float64), 1. +eps(), length=axis_length), range(-eps(Float64), 1. +eps(), length=axis_length)]; kwargs...)
end
using LinearAlgebra
function plot_nullclines!(fig::Figure, p::AbstractWCMNullclineParams, nullcline_axes::AbstractVector{<:AbstractVector};
        xlabel="E", ylabel="I", 
        title="",
        mark_fp=true, arrows_step=nothing
        )
    ax = MakieLayout.Axis(fig, aspect=DataAspect(),
        xlabel=xlabel, ylabel=ylabel,
        xticks=[0,1], yticks=[0,1],
        title=title
    )

    du_defn, dv_defn = field_functions(p)
 
    dus = calculate_field(du_defn, nullcline_axes, p)
    dvs = calculate_field(dv_defn, nullcline_axes, p)
    u_nullclines = Contour.lines(Contour.contour(nullcline_axes..., dus, 0.))
    v_nullclines = Contour.lines(Contour.contour(nullcline_axes..., dvs, 0.))


    if arrows_step !== nothing
        arrow_axis = 0.:arrows_step:1. |> collect
        arrow_axes = [arrow_axis, arrow_axis]
        dus = calculate_field(du_defn, arrow_axes, p)
        dvs = calculate_field(dv_defn, arrow_axes, p)
        strength = vec(sqrt.(dus .^ 2 .+ dvs .^ 2))
        arrows!(ax, arrow_axis, arrow_axis, dus, dvs; arrowcolor=strength, linecolor=strength)
    end

    for line in u_nullclines
        xs, ys = Contour.coordinates(line)
        Makie.lines!(ax, xs, ys; 
            color=:blue
        )
    end

    for line in v_nullclines
        xs, ys = Contour.coordinates(line)
        Makie.lines!(ax, xs, ys, color=:red, 
            linestyle=:dash
        )
    end

    if mark_fp
        fixedpoints = calculate_fixedpoints(p, length(nullcline_axes[1]))
        stability = fixedpoint_stability.(Ref(p), fixedpoints)
        stability_marker = getindex.(Ref(NullclineAnalysis.STABILITY_MARKERS), stability)
        scatter!(ax, Point2f0.(fixedpoints), marker=stability_marker,color=:darkgrey)
    end

    xlims!(ax, 0., 1.)
    ylims!(ax, 0., 1.)

    return ax
end