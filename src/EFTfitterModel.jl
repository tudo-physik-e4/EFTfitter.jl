export EFTfitterModel
export to_correlation_matrix

export get_observables, get_parameters, get_measurements, get_correlations
export get_nuisance_correlations, get_total_covariance, get_covariances, get_measurement_distributions

"""
    struct EFTfitterModel

This is the central type for using EFTfitter.
It comprises all information necessary for performing an analysis.
Only active `Measurement` and `Correlation` objects will be considered.

Fields:  
* `parameters::BAT.NamedTupleDist`
* `measurements::NamedTuple{<:Any, <:Tuple{Vararg{Measurement}}}`
* `measurementdistributions::NamedTuple{<:Any, <:Tuple{Vararg{MeasurementDistribution}}}`
* `correlations::NamedTuple{<:Any, <:Tuple{Vararg{Correlation}}}`
* `nuisances::Union{NamedTuple{<:Any, <:Tuple{Vararg{NuisanceCorrelation}}}, Nothing}`

Constructors:

```julia
EFTfitterModel(
    parameters::BAT.NamedTupleDist,
    measurements::NamedTuple{<:Any, <:Tuple{Vararg{AbstractMeasurement}}},
    correlations::NamedTuple{<:Any, <:Tuple{Vararg{AbstractCorrelation}}},
    nuisances::Union{NamedTuple{<:Any, <:Tuple{Vararg{NuisanceCorrelation}}}, Nothing} = nothing
)
```

Examples:

```julia
model = EFTfitterModel(parameters, measurements, correlations) # no nuisance correlations

model = EFTfitterModel(parameters, measurements, correlations, nuisances) # with nuisance correlations
)
```
"""
struct EFTfitterModel
    parameters::BAT.NamedTupleDist
    measurements::NamedTuple{<:Any, <:Tuple{Vararg{Measurement}}}
    measurementdistributions::NamedTuple{<:Any, <:Tuple{Vararg{MeasurementDistribution}}}
    correlations::NamedTuple{<:Any, <:Tuple{Vararg{Correlation}}}
    nuisances::Union{NamedTuple{<:Any, <:Tuple{Vararg{NuisanceCorrelation}}}, Nothing}
end


function EFTfitterModel(
    parameters::BAT.NamedTupleDist,
    measurements::NamedTuple{<:Any, <:Tuple{Vararg{AbstractMeasurement}}},
    correlations::NamedTuple{<:Any, <:Tuple{Vararg{AbstractCorrelation}}},
    nuisances::Union{NamedTuple{<:Any, <:Tuple{Vararg{NuisanceCorrelation}}}, Nothing} = nothing
)
    measurement_vec, measurement_keys = unpack(measurements)
    correlation_vec, uncertainty_keys = unpack(correlations)

    # convert elements of MeasurementDistribution to Measurement for each bin
    binned_measurements, binned_measurement_keys = convert_to_bins(measurement_vec, measurement_keys)
    # use only active measurements/bins
    active_measurements, active_measurement_keys, corrs = only_active_measurements(binned_measurements, binned_measurement_keys, correlation_vec)
    # use only active uncertainties and correlations
    active_measurements, active_correlations, uncertainty_keys = only_active_uncertainties(active_measurements, corrs, uncertainty_keys)

    correlation_nt = namedtuple(uncertainty_keys, active_correlations)
    measurement_nt = namedtuple(active_measurement_keys, active_measurements)
    meas_dists_nt  = create_distributions(measurements, uncertainty_keys)
    nuisances_nt   = only_active_nuisances(nuisances, active_measurement_keys, uncertainty_keys)
    
    params = add_nuisance_parameters(parameters, nuisances_nt)

    return EFTfitterModel(params, measurement_nt, meas_dists_nt, correlation_nt, nuisances_nt)
end
#=============================================================#
"""
    get_parameters(m::EFTfitterModel)

Returns model parameters.
"""
get_parameters(m::EFTfitterModel) = m.parameters


"""
    get_measurements(m::EFTfitterModel)

Returns a `NamedTuple` with the `Measurement`s in the `EFTfitterModel`.
"""
get_measurements(m::EFTfitterModel) = m.measurements


"""
    get_measurement_distributions(m::EFTfitterModel)

Returns a `NamedTuple` with the `MeasurementDistribution`s in the `EFTfitterModel`.
"""
get_measurement_distributions(m::EFTfitterModel) = m.measurementdistributions


"""
    get_correlations(m::EFTfitterModel)

Returns a `NamedTuple` with the `Correlations`s in the `EFTfitterModel`.
"""
get_correlations(m::EFTfitterModel) = m.correlations


"""
    get_nuisance_correlations(m::EFTfitterModel)

Returns a `NamedTuple` with the `NuisanceCorrelation`s in the `EFTfitterModel`.
"""
get_nuisance_correlations(m::EFTfitterModel) = m.nuisances


"""
    get_observables(m::EFTfitterModel)

Returns a `NamedTuple` with the `Observable`s in the `EFTfitterModel`.
Note: The upper and lower limits are ignored and for each unique `Function`s only one `Observable` is returned.
"""
function  get_observables(model::EFTfitterModel)
    meas = get_measurements(model)
    obs = unique(Observable.([m.observable.func for m in values(meas)]))
    obs_names = [string(o.func) for o in obs]
    observables_nt = namedtuple(obs_names, obs)
end



#========================================================#

function unpack(nt::NamedTuple{N,T}) where {N,T}
     return collect(nt), collect(keys(nt))
end


function convert_to_bins(
    measurements::Array{<:AbstractMeasurement, 1},
    measurement_keys::Array{Symbol, 1}
)
    binned_measurements = reduce(vcat, convert_to_bins.(measurements))
    binned_meas_keys = reduce(vcat, keys_of_bins.(measurements, measurement_keys))

    if isa(binned_measurements, Measurement) # when only one measurement is active
        return [binned_measurements], [binned_meas_keys]
    end

    return binned_measurements, binned_meas_keys
end


function convert_to_bins(m::Measurement)
    return m
end


function convert_to_bins(md::MeasurementDistribution)
    nbins = length(md.value)
    uncertainties = [[u[i] for u in md.uncertainties] for i in 1:nbins]

    unc_nt = namedtuple.(Ref(keys(md.uncertainties)), uncertainties)
    return Measurement.(md.observable, md.value, unc_nt, md.active)
end



function only_active_measurements(
    measurements::Array{Measurement, 1},
    measurement_keys::Array{Symbol, 1},
    correlations::Array{<:AbstractCorrelation, 1}
)
    active_idxs = filter(x->measurements[x].active, 1:length(measurements))
    measurements = measurements[active_idxs]
    measurement_keys = measurement_keys[active_idxs]

    corrs = get_correlations.(correlations, Ref(active_idxs))

    return measurements, measurement_keys, corrs
end

#TODO: rename?
function get_correlations(
    correlation::Correlation,
    idxs::Array{<:Integer, 1}
)
    n=length(idxs)
    # ignore size of given correlation matrix if active==false
    if correlation.active == false
        return Correlation(zeros(n, n), false)
    end

    return Correlation(correlation.matrix[idxs, idxs], correlation.active)
end

function get_correlations(
    correlation::NoCorrelation,
    idxs::Array{<:Integer, 1}
)
    n = length(idxs)
    return Correlation(Matrix{Float64}(I, n, n), correlation.active)
end



#--------------------
function create_distributions(
    m::NamedTuple{<:Any, <:Tuple{Vararg{AbstractMeasurement}}},
    uncertainty_keys::Vector{Symbol}
)
    active_idxs = [i for i in 1:length(m) if (isa(m[i], MeasurementDistribution) && any(m[i].active))]
    dists = [only_active_bins(m[i], uncertainty_keys) for i in active_idxs]
    dists_keys = [keys(m)[i] for i in active_idxs]

    if isempty(dists)
        return NamedTuple()
    else
        return namedtuple(dists_keys, dists)
    end
end

function only_active_bins(md::MeasurementDistribution, uncertainty_keys::Vector{Symbol})
    obs = md.observable[md.active]
    vals = md.value[md.active]

    unc = [md.uncertainties[u][md.active] for u in uncertainty_keys]
    uncertainties = namedtuple(uncertainty_keys, unc)
    bin_names = md.bin_names[md.active]

    return MeasurementDistribution(obs, vals, uncertainties=uncertainties, bin_names=bin_names)
end


function keys_of_bins(m::Measurement, key::Symbol)
    return key
end

function keys_of_bins(md::MeasurementDistribution, key::Symbol)
    return [Symbol(String(key)*"_"*String(b)) for b in md.bin_names]
end



function only_active_uncertainties(
    measurements::Array{Measurement, 1},
    correlations::Array{<:AbstractCorrelation, 1},
    uncertainty_keys::Array{Symbol, 1}
)
    active_idxs = filter(x->correlations[x].active, 1:length(correlations))
    active_correlations = correlations[active_idxs]
    uncertainty_keys = uncertainty_keys[active_idxs]

    measurements = [@set m.uncertainties = select(m.uncertainties, uncertainty_keys) for m in measurements]

    return measurements, active_correlations, uncertainty_keys
end


function only_active_nuisances(
    nuisances::Nothing,
    active_meas_keys::Array{Symbol, 1},
    unc_keys::Array{Symbol, 1}
)
    return nothing
end


function only_active_nuisances(
    nuisances::NamedTuple{<:Any, <:Tuple{Vararg{NuisanceCorrelation}}},
    active_meas_keys::Array{Symbol, 1},
    unc_keys::Array{Symbol, 1}
)
    nuisance_vec, nuisance_keys = unpack(nuisances)

    active_nuisances = NuisanceCorrelation[]
    active_nuisance_keys = Symbol[]

    for (nui, nui_k) in zip(nuisance_vec, nuisance_keys)
        if check_nuisance(nui, active_meas_keys, unc_keys)
            push!(active_nuisances, nui)
            push!(active_nuisance_keys, nui_k)
        end
    end

    if length(active_nuisances) > 0
        return namedtuple(active_nuisance_keys, active_nuisances)
    else
        return nothing
    end
end


function check_nuisance(
    nuisance::NuisanceCorrelation,
    active_meas_keys::Array{Symbol, 1},
    unc_keys::Array{Symbol, 1}
)
    if nuisance.meas1 == nuisance.meas2
        return false
    elseif !(nuisance.unc_key in unc_keys)
        return false
    elseif !(nuisance.meas1 in active_meas_keys)
        return false
    elseif !(nuisance.meas2 in active_meas_keys)
        return false
    end

    return true
end



function get_total_covariance(m::EFTfitterModel)
    covs = get_covariances(m)
    total_cov = sum(covs)
    total_cov = convert(Array{Float64,2}, total_cov)

    return total_cov
end


function get_covariances(m::EFTfitterModel)
    unc_values = [[meas.uncertainties[u] for meas in m.measurements] for u in keys(m.correlations)]
    corrs = [c.matrix for c in m.correlations]

    covs = [σ*ρ*σ for (σ, ρ) in zip(diagm.(unc_values), corrs)]
end



function Nbins(measurements::NamedTuple{<:Any, <:Tuple{Vararg{AbstractMeasurement}}})
    sum([length(m.value) for m in measurements])
end



function get_indices(measurements, symb)
    meas_vec, meas_keys = unpack(measurements)
    names = keys_of_bins(meas_vec, meas_keys)
    ind = [findfirst(x -> x==symb, names)]
    
    if ind[1] != nothing
        return ind
    else
        ind = [findfirst(x -> x==Symbol(String(symb)*"_bin$i"), names) for i in 1:length(measurements[symb].value)]
    end
end



"""
    to_correlation_matrix(
        measurements::NamedTuple{<:Any, <:Tuple{Vararg{AbstractMeasurement}}},
        correlations::Tuple{Symbol,Symbol, Union{<:Real, Array{<:Real, 2}}}...
    )

    Returns a `Matrix{Float64}` as the correlation matrix for the measurements with all diagonal-elements being unity.
    The correlations are specified by passing tuples of two `Symbols` (the keys of the `Measurements`) with a value or matrix for the correlation coefficients.
    If the matrix is not positive-definite, a warning is shown but the matrix is still returned.
    
    Example:
    ```julia
    dist_corr = [1.0 0.5 0.0;
                 0.5 1.0 0.0;
                 0.0 0.0 1.0]
    
    another_corr_matrix = to_correlation_matrix(
        measurements,
        (:Meas1, :Meas2, 0.4), 
        (:Meas1, :MeasDist, 0.1), 
        (:MeasDist, :MeasDist, dist_corr), 
        (:MeasDist_bin2, :MeasDist_bin3, 0.3),
    )
    ```
"""
function to_correlation_matrix(
    measurements::NamedTuple{<:Any, <:Tuple{Vararg{AbstractMeasurement}}},
    correlations::Any...#Tuple{Symbol,Symbol, Union{<:Real, Array{<:Real, 2}}}... #TODO: typing
)
    nmeas = Nbins(measurements)
    corr = Matrix{Float64}(I, nmeas, nmeas)

    for c in correlations #TODO: split into functions
        if isa(c, Tuple{Symbol,Symbol, Union{<:Real, Array{<:Real, 2}}})
            i = get_indices(measurements, c[1])
            j = get_indices(measurements, c[2])

            corr[i, j] .= c[3]
            corr[j, i] .= c[3]
            
        elseif isa(c, Tuple{Any, Union{<:Real, Array{<:Real, 2}}})
            i = reduce(vcat, get_indices.(Ref(measurements), c[1]))
            corr[i, i] .= c[2]
        end
    end

    corr[diagind(corr)] .= 1

    if !isposdef(corr)
        @warn "The correlation matrix $corr is not positive semidefinite"
    end

    return corr
end



function has_nuisance_correlations(m::EFTfitterModel)
    m.nuisances == nothing ? (return false) : (return true)
end


function add_nuisance_parameters(
    parameters::NamedTupleDist, 
    nuisances_nt::NamedTuple{<:Any, <:Tuple{Vararg{NuisanceCorrelation}}}
)
    param_keys = [keys(parameters)...,  keys(nuisances_nt)...]
    
    nuisances = [nui.prior for nui in nuisances_nt]
    param_dists = [values(parameters)...,  nuisances...]
    nt = namedtuple(param_keys, param_dists)
    
    return BAT.NamedTupleDist(nt)
end

function add_nuisance_parameters(
    parameters::NamedTupleDist, 
    nuisances_nt::Nothing
)
    return parameters
end



# Todo: rename?
function keys_of_bins(
    measurements::Array{T, 1},
    measurement_keys::Array{Symbol, 1}
) where T <: AbstractMeasurement
    keys = reduce(vcat, keys_of_bins.(measurements, measurement_keys))
    return keys
end
