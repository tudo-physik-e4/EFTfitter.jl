# plot recipe for plotting a Measurement object
@recipe function f(
    measurement::Measurement;
    uncertainties = true
)
    total_unc = if uncertainties == false
        0.0
    elseif uncertainties == true
        total_uncertainty(measurement, keys(measurement.uncertainties))
    else
        total_uncertainty(measurement, uncertainties)
    end

    @series begin
        seriestype --> :hline
        ribbon --> total_unc
        [measurement.value]
    end

end

# plot recipe for plotting a MeasurementDistribution object
@recipe function f(
    measurement::MeasurementDistribution;
    uncertainties = true,
    bin_names = String.(measurement.bin_names)
)
    n = length(measurement.value)
    xs = 1:n+1
    ys = [measurement.value[1], measurement.value...]

    total_uncs = if uncertainties == false
        zeros(n)
    elseif uncertainties == true
        total_uncertainty(measurement, keys(measurement.uncertainties))
    else
        total_uncertainty(measurement, uncertainties)
    end

    seriestype = get(plotattributes, :seriestype, :steppre)

    #TODO adapt for arbitrary bin width
    if seriestype == :steppre
        @series begin
            seriestype := :steppre
            ribbon --> [1, total_uncs...]
            label --> "measurements"
            xticks -->  (collect(xs).+0.5, bin_names)
            xs, ys
        end

    elseif seriestype == :scatter
        @series begin
            seriestype := :scatter
            xerror --> 0.5
            yerror --> total_uncs
            seriescolor --> :orange
            label --> "measurements"
            xticks --> (collect(xs).+0.5, bin_names)
            xs[2:end].-0.5, ys[2:end]
        end
    end
end

# plot recipe for plotting a particular Measurement from a NamedTuple
@recipe function f(
    measurements::NamedTuple{<:Any, <:Tuple{Vararg{AbstractMeasurement}}},
    name::Symbol;
    uncertainties = true
)
    @series begin
        label --> string(name)
        uncertainties --> uncertainties
        measurements[name]
    end
end

# plot recipe for plotting a particular MeasurementDistribution from a NamedTuple
@recipe function f(
    measurements::NamedTuple{<:Any, <:Tuple{Vararg{AbstractMeasurement}}},
    name::Symbol;
    uncertainties = true
)
    @series begin
        label --> string(name)
        uncertainties --> uncertainties
        measurements[name]
    end
end


function total_uncertainty(
    measurement::Union{Measurement, MeasurementDistribution},
    uncertainty::Symbol
)
    return measurement.uncertainties[uncertainty]
end

# calculate the total uncertainty of a measurement by adding the quadratures
function total_uncertainty(
    measurement::Union{Measurement, MeasurementDistribution},
    uncertainty::Union{Tuple{Vararg{Symbol}}, Vector{Symbol}}
)
    total = zeros(length(measurement.uncertainties[1]))
    for u in uncertainty
        total .+= measurement.uncertainties[u].^2
    end
    return sqrt.(total)
end
