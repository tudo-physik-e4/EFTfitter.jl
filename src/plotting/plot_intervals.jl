# strings of parameter names that are not fixed
function get_parameter_names(samples::DensitySampleVector)
    parameter_names = collect(string.(keys(BAT.varshape(samples))))
    parameter_names = [parameter_names[i] for i in 1:length(parameter_names) if BAT.varshape(samples)._accessors[i].len != 0]
    return parameter_names
end


@recipe function f(
    maybe_shaped_samples::DensitySampleVector,
    p::Float64;
    parameter_names = get_parameter_names(maybe_shaped_samples),
    y_positons = collect(1:length(parameter_names))*-1,
    y_offset = 0,
    bins = 200,
    atol = 0,
)
    samples = BAT.unshaped.(maybe_shaped_samples)
    
    for i in 1:length(y_positons)

        lower, upper =  get_interval_edges(samples, i, p, bins=bins, atol=atol)
        w = 0.5*(upper-lower)
        middle = lower .+ w
        
        ypos = fill(y_positons[i], length(middle)).-y_offset
        
        if i == 1
            intern_label = get(plotattributes, :label, "")
        else
            intern_label = ""
        end
        
        internal_color = get(plotattributes, :markerstrokecolor , :steelblue)
        
        @series begin
            seriestype --> :path
            markerstrokealpha := 0
            linewidth := 0
            markerstrokewidth --> 1
            markersize --> 5
            markerstrokecolor --> :steelblue
            linecolor := internal_color
            xerror := w
            yticks --> (y_positons, parameter_names)
            label := intern_label
            xguide --> "Ci"

            (middle, ypos)
        end
    end
end
