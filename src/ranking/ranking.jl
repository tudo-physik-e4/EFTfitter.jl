include("ranking_criteria.jl")

abstract type AbstractRanks end

"""
    struct MeasurementRanks

Type for the results of the measurement ranking.

Fields:  
* `names::Vector{Symbol}`: Names of the measurements.
* `values::Vector{Float64}`: Values of the ranking, calculated according to the criterion.
* `criterion<:AbstractRankingCriterion`: Criterion that was used for ranking.

Constructors:

```julia
function MeasurementRanks(
    names::Vector{Symbol}, 
    values::Vector{Float64}, 
    criterion::AbstractRankingCriterion;
    order::Symbol = :values,
    rev::Bool = false
)
```
Keyword arguments:
* `order::Symbol = :values`: Specifies how the `names` and `values` are sorted. Further options are 
`:names` for alphabetical sorting based on the measurement names or `:none` for keeping the initial order of measurements. 
* `rev::Bool=false`: Switch to invert the order.
"""
struct MeasurementRanks <: AbstractRanks
    names::Vector{Symbol}
    values::Vector{Float64}
    criterion::AbstractRankingCriterion
    
    function MeasurementRanks(
        names::Vector{Symbol}, 
        values::Vector{Float64}, 
        criterion::AbstractRankingCriterion;
        order::Symbol = :values,
        rev::Bool = false
    )
        sorted_names, sorted_values = sort_ranking_results(order, names, values, rev=rev)
        return new(sorted_names, sorted_values, criterion)
    end
end


"""
    struct UncertaintyRanks

Type for the results of the unceertainty ranking.

Fields:  
* `names::Vector{Symbol}`: Names of the uncertainty categories.
* `values::Vector{Float64}`: Values of the ranking, calculated according to the criterion.
* `criterion<:AbstractRankingCriterion`: Criterion that was used for ranking.

Constructors:

```julia
function UncertaintyRanks(
    names::Vector{Symbol}, 
    values::Vector{Float64}, 
    criterion::AbstractRankingCriterion;
    order::Symbol = :values,
    rev::Bool = false
)
```
Keyword arguments:
* `order::Symbol = :values`: Specifies how the `names` and `values` are sorted. Further options are 
`:names` for alphabetical sorting based on the uncertainty names or `:none` for keeping the initial order of uncertainty categories. 
* `rev::Bool=false`: Switch to invert the order.
"""
struct UncertaintyRanks <: AbstractRanks
    names::Vector{Symbol}
    values::Vector{Float64}
    criterion::AbstractRankingCriterion
    
    function UncertaintyRanks(
        names::Vector{Symbol}, 
        values::Vector{Float64}, 
        criterion::AbstractRankingCriterion;
        order::Symbol = :values,
        rev::Bool = false
    )
        sorted_names, sorted_values = sort_ranking_results(order, names, values, rev=rev)
        return new(sorted_names, sorted_values, criterion)
    end
end



"""
    rank_measurements(
        model::EFTfitterModel;
        sampling_algorithm::BAT.AbstractSamplingAlgorithm = MCMCSampling(),
        criterion = SmallestIntervalsSum(p=0.9),
        order = :values,
        rev = true
    )
    
Computes a ranking of the individual measurement in the `EFTfitterModel` based
on a certain ranking criterion by performing the sampling according to the specified `sampling_algorithm`.
Returns a `MeasurementRank` object.

By default, the summed width of all one-dimensional marginalized smallest intervals containing 90% of the posterior
probability is used as the ranking criterion. 
The measurements are ranked by the relative increase of this value.

Available ranking criteria:

* `SmallestIntervalsSum(p=0.9)`: summed width of all one-dimensional marginalized smallest intervals containing p=90% of the posterior probability
* `SmallestInterval(key=:C1, p=0.9)`: width of the one-dimensional marginalized smallest interval of parameter `:C1` containing p=90% of the posterior probability.
"""
function rank_measurements(
    model::EFTfitterModel;
    sampling_algorithm::BAT.AbstractSamplingAlgorithm = MCMCSampling(),
    criterion = SumOfSmallestIntervals(),
    order = :values,
    rev = true
)
    @info "Begin ranking of measurements"
    original = criterion_value(criterion, model, sampling_algorithm=sampling_algorithm)

    deactivated_models, measurement_keys = measurement_models(model)
    values = [criterion_value(criterion, m, sampling_algorithm=sampling_algorithm) for m in deactivated_models]

    relative_increase = calculate_ranks(original, values)

    return MeasurementRanks(collect(measurement_keys), relative_increase, criterion, order=order, rev=rev)
end


"""
    rank_uncertainties(
        model::EFTfitterModel;
        sampling_algorithm::BAT.AbstractSamplingAlgorithm = MCMCSampling(),
        criterion = SmallestIntervalsSum(p=0.9),
        order = :values,
        rev = true
    )
    
Computes a ranking of the types of uncertainty in the `EFTfitterModel` based
on a certain ranking criterion by performing the sampling according to the specified `sampling_algorithm`.
Returns a `MeasurementRank` object.

By default, the summed width of all one-dimensional marginalized smallest intervals containing 90% of the posterior
probability is used as the ranking criterion. 
The uncertainty types are ranked by the relative decrease of this value.

Available ranking criteria:

* `SmallestIntervalsSum(p=0.9)`: summed width of all one-dimensional marginalized smallest intervals containing p=90% of the posterior probability
* `SmallestInterval(key=:C1, p=0.9)`: width of the one-dimensional marginalized smallest interval of parameter `key=:C1` containing p=90% of the posterior probability
"""
function rank_uncertainties(
    model::EFTfitterModel;
    sampling_algorithm::BAT.AbstractSamplingAlgorithm = MCMCSampling(),
    criterion = SumOfSmallestIntervals(),
    order = :values,
    rev = true
)
    @info "Begin ranking of uncertainties"
    original = criterion_value(criterion, model, sampling_algorithm=sampling_algorithm)

    deactivated_models, uncertainty_keys = uncertainty_models(model)
    values = [criterion_value(criterion, m, sampling_algorithm=sampling_algorithm) for m in deactivated_models]
    relative_decrease = -1*calculate_ranks(original, values)

    return UncertaintyRanks(collect(uncertainty_keys), relative_decrease, criterion, order=order, rev=rev)
end


"""
    function measurement_models(model::EFTfitterModel)
    
Creates a `Vector` of `EFTfitterModel` where always one of the initially active measurements is deactivated at a time.
Returns a `Vector` of `EFTfitterModel` and a `Vector` of `Symbol` with the names of the currently deactivated measurement.
"""
function measurement_models(model::EFTfitterModel)
    m = model.measurements
    measurement_keys = collect(keys(m))

    if length(measurement_keys) < 2
         throw(ArgumentError("Only one measurement active. Ranking of measurements not possible."))
    end

    models = Vector{EFTfitterModel}(undef, length(m))
    for i in 1:length(m)
        meas = collect(values(model.measurements))
        meas[i] = deactivate_measurement(meas[i])
        meas_nt = namedtuple(measurement_keys, meas)
        models[i] = EFTfitterModel(model.parameters, meas_nt, model.correlations)
    end

    return models, measurement_keys
end


"""
    function uncertainty_models(model::EFTfitterModel)
    
Creates a `Vector` of `EFTfitterModel` where always one of the initially active uncertainty types is deactivated at a time.
Returns a `Vector` of `EFTfitterModel` and a `Vector` of `Symbol` with the names of the currently deactivated uncertainty types.
"""
function uncertainty_models(model::EFTfitterModel)
    cors = model.correlations
    uncertainty_keys = collect(keys(cors))

    if length(uncertainty_keys) < 2
         throw(ArgumentError("Only one type of uncertainties active. Ranking of uncertainties not possible."))
    end

    models = Vector{EFTfitterModel}(undef, length(cors))
    for i in 1:length(cors)
        cors = collect(values(model.correlations))
        cors[i] = deactivate_uncertainty(cors[i])
        cors_nt = namedtuple(uncertainty_keys, cors)
        models[i] = EFTfitterModel(model.parameters, model.measurements, cors_nt)
    end

    return models, uncertainty_keys
end


function deactivate_measurement(meas::Measurement)
    return Measurement(meas.observable, meas.value, meas.uncertainties, false)
end

function deactivate_uncertainty(cor::Correlation)
    return Correlation(cor.matrix, false)
end


"""
    function criterion_value(
        criterion::AbstractRankingCriterion,
        model::EFTfitterModel;
        sampling_algorithm::BAT.AbstractSamplingAlgorithm=MCMCSampling()
    )
    
Computes the value of the `criterion` for the `model` using the `sampling_algorithm`.

Available ranking criteria:

* `SmallestIntervalsSum(p=0.9)`: summed width of all one-dimensional marginalized smallest intervals containing p=90% of the posterior probability
* `SmallestInterval(key=:C1, p=0.9)`: width of the one-dimensional marginalized smallest interval of parameter `key=:C1` containing p=90% of the posterior probability
"""
function criterion_value(
    criterion::AbstractRankingCriterion,
    model::EFTfitterModel;
    sampling_algorithm::BAT.AbstractSamplingAlgorithm=MCMCSampling()
)
    posterior = PosteriorDensity(EFTfitterDensity(model), model.parameters)
    samples = BAT.bat_sample(posterior, sampling_algorithm).result
    criterion_value = apply_criterion(criterion, samples)
end


function calculate_ranks(original_value::Real, values::AbstractArray{<:Real})
    return relative_increase = (values.-original_value)./original_value
end

function calculate_ranks(original_values::AbstractArray{<:Real}, values::AbstractArray{<:AbstractArray{<:Real, 1}, 1})
    relative_increase = [sum((val.-original_values)./original_values) for val in values]
end


function sort_ranking_results(
    order::Symbol, 
    names::Vector{Symbol}, 
    values::Vector{Float64};
    rev=false
)
    if order == :names
        sorted = sortperm(names, rev=rev)
        return names[sorted], values[sorted]
        
    elseif order == :values
        sorted = sortperm(values, rev=rev)
        return names[sorted], values[sorted]
        
    elseif order == :none
        return names, values

    else 
        return names, values
    end
end
