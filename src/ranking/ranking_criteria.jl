export SumOfSmallestIntervals
export SmallestInterval
export HighestDensityRegion

export get_smallest_interval_edges

abstract type AbstractRankingCriterion end

"""
    struct SmallestInterval <: AbstractRankingCriterion

Type for specifying the ranking criterion to be the total width of the
one-dimensional marginalized smallest interval containing a probability mass `p` for a certain parameter.

Fields:  
* `key::Symbol`: Name of the parameter that is used for the marginalized intervals.
* `p::Float64 = 0.9`: Probability mass to be enclosed in the smallest interval.

Constructors:

```julia
SmallestInterval(key::Symbol, p=0.9, bins=200)
```
"""
@with_kw struct SmallestInterval{T} <: AbstractRankingCriterion
    key::Symbol
    p::Float64 = 0.9
    bins::T = 200
end


"""
    struct SumOfSmallestIntervals <: AbstractRankingCriterion

Type for specifying the ranking criterion to be the summed width of all 
one-dimensional marginalized smallest intervals containing a probability mass `p`.

Fields:  
* `p::Float64 = 0.9`: Probability mass to be enclosed in the smallest intervals.
* `bins::T = 200`: Number of bins for the histograms to calculate interval widths. 

Constructors:

```julia
SumOfSmallestIntervals(p=0.9, bins=200)
```
"""
@with_kw struct SumOfSmallestIntervals{T} <: AbstractRankingCriterion
    p::Float64 = 0.9
    bins::T = 200
end



"""
    struct HighestDensityRegion <: AbstractRankingCriterion

Type for specifying the ranking criterion to be the area of the two-dimensional 
marginal posterior region containing a probability mass `p` for two certain parameters.

Fields:  
* `keys::NTuple{2, Symbol}`: Names of the parameters that are used for the marginalized regions.
* `p::Float64 = 0.9`: Probability mass to be enclosed in the smallest interval.
* `bins::T = 200`: Number of bins for the histograms to calculate interval widths. 

Constructors:

```julia
HighestDensityRegion(keys::NTuple{2, Symbol}, p=0.9, bins=200)
```
"""
@with_kw struct HighestDensityRegion{N, T} <: AbstractRankingCriterion
    keys::NTuple{N, Symbol}
    p::Float64 = 0.9
    bins::T = 200
end

"""
    apply_criterion(criterion<:AbstractRankingCriterion, samples::DensitySampleVector)
    
Applies the `criterion` to the `samples` and returns the corresponding value.

Available ranking criteria:

* `SumOfSmallestIntervals(p=0.9)`: summed width of all one-dimensional marginalized smallest intervals containing p=90% of the posterior probability
* `SmallestInterval(key=:C1, p=0.9)`: width of the one-dimensional marginalized smallest interval of parameter `:C1` containing p=90% of the posterior probability.
* `HighestDensityRegion(keys::NTuple{2, Symbol}, p=0.9, bins=200): volume of the highest density regione
"""
function apply_criterion(criterion::SumOfSmallestIntervals, samples::DensitySampleVector)
    return get_all_intervals(samples, criterion.p, bins=criterion.bins)
end

function apply_criterion(criterion::SmallestInterval, samples::DensitySampleVector)
    return get_volume(samples, (criterion.key, ), criterion.p, bins=criterion.bins)
end

function apply_criterion(criterion::HighestDensityRegion, samples::DensitySampleVector)
    return get_volume(samples, criterion.keys, criterion.p, bins=criterion.bins)
end

#------------------------------------------------------------------
# helper functions for calculating the intervals 

# total width of 1d marginal smallest intervals
function get_interval(samples, key, p; bins=200)
    marg = BAT.MarginalDist(samples, key, bins=bins)

    marghist = convert(Histogram, marg.dist)
    hist_p = BAT.get_smallest_intervals(marghist, [p])[1][1]

    lower, upper = BAT.get_smallest_interval_edges(hist_p)
    return sum(upper .- lower)
end



"""
    get_smallest_interval_edges(
        samples::DensitySampleVector, 
        key::Union{Symbol, Real}, 
        p::Real; 
        bins=200,
        atol=0.0)
    
Calculates the edges of the smallest intervals containing the fraction `p` of the probability of the marginal distribution of parameter `key`.

Returns a `NamedTuple` with the keys `lower` and `upper` which both contain `Vector`s with the corresponding bin edges.
Keywords:
* `bins=200`: The number of bins used fpr calculating the intervals and edges.
* `atol=0.0`: Intervals are joined together when they are seperated by less than this value.
"""
function get_smallest_interval_edges(
    samples::DensitySampleVector, 
    key::Union{Symbol, Real}, 
    p::Real; 
    bins=200,
    atol=0.0)
    
    marg = BAT.MarginalDist(samples, key, bins=bins)

    marghist = convert(Histogram, marg.dist)
    hist_p = BAT.get_smallest_intervals(marghist, [p])[1][1]

    lower, upper = BAT.get_interval_edges(hist_p, atol=atol)
    return (lower=lower, upper=upper)
end


function get_volume(samples, keys, p; bins=200)
    marghist = get_marginal_hist(samples, keys, bins=bins)
    
    hist_p = BAT.get_smallest_intervals(marghist, [p])[1][1]
    
    vol = bin_volume(hist_p)
    n = count(x->x>0, hist_p.weights) 
    return vol*n
end

function bin_volume(hist::StatsBase.Histogram)
    edges = hist.edges
    vol = prod([abs(edges[i][1]-edges[i][2]) for i in 1:length(edges)])
    return vol
end


function get_all_intervals(samples, p; bins=200)
    sample_keys = BAT.all_active_names(BAT.varshape(samples))

    intervals = [get_volume(samples, (Meta.parse(k),), p, bins=bins) for k in sample_keys]
    return intervals
end


function get_marginal_hist(
    maybe_shaped_samples::BAT.DensitySampleVector,
    key::Union{NTuple{n,Integer}, NTuple{n,Union{Symbol, Expr}}} where n;
    bins = 200,
    closed::Symbol = :left,
    filter::Bool = false
)
    samples = unshaped.(maybe_shaped_samples)

    if filter
        samples = BAT.drop_low_weight_samples(samples)
    end

    idxs = BAT.asindex.(Ref(maybe_shaped_samples), key)
    s = Tuple(BAT.flatview(samples.v)[i, :] for i in idxs)

    edges = if isa(bins, Integer)
        BAT._get_edges(s, (bins,), closed)
    else
        Tuple(BAT._get_edges(s[i], bins[i], closed) for i in 1:length(bins))
    end

    return fit(Histogram, s, FrequencyWeights(samples.weight), edges, closed = closed)
end
