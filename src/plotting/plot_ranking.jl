@recipe function f(
    result::AbstractRanks
)
    names = string.(result.names)
    ylabel = rankingcriterion(result)
    
    # prevent pyplot from converting underscores to subscripts
    if _plots_module().backend() == _plots_module().PyPlotBackend()
        names = replace.(names, "_" => "-")
        ylabel = replace(ylabel, "_" => "-")
    end
    
    @series begin
        seriestype --> :bar
        legend --> false
        yguide --> ylabel
        xrotation --> 0

        names, result.values
    end
end



function rankingcriterion(result::MeasurementRanks)
    if isa(result.criterion, SumOfSmallestIntervals)
        return "sum of rel. increases of smallest "*string(result.criterion.p*100)*"% intervals"
        
    elseif isa(result.criterion, SmallestInterval)
        return "rel. increase of smallest "*string(result.criterion.p*100)*"% interval of "*string(result.criterion.key)
    
    elseif isa(result.criterion, HighestDensityRegion)
        return "rel. increase of "*string(result.criterion.p*100)*"% HDR("*join(map(string, result.criterion.keys), ", ")*")"
    end
        
    return ""
end


function rankingcriterion(result::UncertaintyRanks)
    if isa(result.criterion, SumOfSmallestIntervals)
        return "sum of rel. decreases of smallest "*string(result.criterion.p*100)*"% intervals"
        
    elseif isa(result.criterion, SmallestInterval)
        return "rel. decrease of smallest "*string(result.criterion.p*100)*"% interval of "*string(result.criterion.key)
    
    elseif isa(result.criterion, HighestDensityRegion)
        return "rel. decrease of "*string(result.criterion.p*100)*"% HDR of "*string(result.criterion.keys)
    end

    return ""
end
