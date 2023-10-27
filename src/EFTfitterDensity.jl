export EFTfitterDensity

#------------- Nuisance Correlation Indices -----------------------------------------#
# internal struct for nuisance correlations (not to be confused with user-interface NuisanceCorrelation)
# TODO: rename ?
struct NuisanceCorrelationIndices
    unc::Int # index of uncertainty 
    i::Int # index of measurement 1
    j::Int # index of measurement 2
    key::Symbol # parameter name for this nuisance correlation
end


#------------- Covariance Matrix Type Indicator -----------------------------------------#
# internal type for dispatch  
abstract type CovOrInverseMatrix end

# when covariance matrix is changed during sampling,
# i.e. when nuisance correlations or model uncertainties need to be added
struct CovMatrix{T} <: CovOrInverseMatrix 
    m::T
    diagonal::Vector{Float64}
end

# when covariance matrix is fix and not changed during sampling,
# i.e. when no nuisance correlations or model uncertainties are used
struct InvCovMatrix{T} <: CovOrInverseMatrix 
    m::T
    diagonal::Vector{Float64}
end


#------------- Model Uncertainties Type Indicator ---------------------------------------#
abstract type ModelUncertaintiesIndicator end

struct HasModelUncertainties <: ModelUncertaintiesIndicator end
struct NoModelUncertainties <: ModelUncertaintiesIndicator end

get_model_uncertainties_indicator(predictions::Vector{<:Real}) = NoModelUncertainties()
get_model_uncertainties_indicator(predictions) = HasModelUncertainties()


#------------- Nuisance Correlations Indicator -----------------------###----------------#
abstract type NuisanceCorrelationsIndicator end
struct NoNuissanceCorrelations <: NuisanceCorrelationsIndicator end

struct NuisanceCorrelations{T} <: NuisanceCorrelationsIndicator 
    nuisance_correlations::Vector{NuisanceCorrelationIndices}
    covs::Vector{T}
end

NuisanceCorrelations(nuisances, m::EFTfitterModel) = NuisanceCorrelations(nuisances, get_covariances(m))


#------------- Limits Indicator ---------------------------------------------------------#
# Status indicating whether upper limits are used or not
abstract type LimitsIndicator end
struct HasLimits <: LimitsIndicator end
struct NoLimits <: LimitsIndicator end
# TODO: Limits not yet used

#------------- Bounds Check Indicator ---------------------------------------------------------#
# Status indicating whether observable bounds should be checked or not
abstract type BoundsCheckIndicator end
struct BoundsCheck <: BoundsCheckIndicator end
struct NoBoundsCheck <: BoundsCheckIndicator end



#------------- EFTfitterDensity ---------------------------------------------------------#
struct EFTfitterDensity{
    M<:CovOrInverseMatrix,
    B<:BoundsCheckIndicator, 
    MU<:ModelUncertaintiesIndicator, 
    NC<:NuisanceCorrelationsIndicator, 
    L<:LimitsIndicator
}
    measured_values::Vector{Float64}
    observable_functions::Vector{Function}
    observable_mins::Vector{Float64}
    observable_maxs::Vector{Float64}
    observable_weights::Vector{Float64} # observable weights 
    crossmatrix::M # CovMatrix or InvCovMatrix
    check_bounds::B # whether to check observable bounds or not
    predictions::Matrix{Float64} 
    prediction_uncertainties::Matrix{Float64} # only used if model uncertainties are present
    # limit_distributions::Vector{Distribution} # only used if limits are present
    # limit_functions::Vector{Function}      # only used if limits are present
    # limit_predictions::Matrix{Float64}    # only used if limits are present
    # limit_uncertainties::Matrix{Float64}  # only used if limits are present
    model_uncertainties::MU  # HasModelUncertainties or NoModelUncertainties
    limits::L   # HasLimits or NoLimits
    nuisance_correlations::NC   # NuisanceCorrelations or NoNuissanceCorrelations
end

@inline DensityInterface.DensityKind(::EFTfitterDensity) = IsDensity()


"""
    get_crossmatrix(mus::NoModelUncertainties, ncs::NoNuissanceCorrelations, m::EFTfitterModel, weights)

For models without model uncertainties and without nuisance correlations, this function 
computes the cross matrix using the observable weights included in the matrix.
"""
function get_crossmatrix(mus::NoModelUncertainties, ncs::NoNuissanceCorrelations, m::EFTfitterModel, weights)
    invcov_weighted = weights .* inv(get_total_covariance(m))
    M_invcov = m.CovarianceType(invcov_weighted)

    return InvCovMatrix(M_invcov, diag(invcov_weighted))
end


"""
    get_crossmatrix(mus::ModelUncertaintiesIndicator, ncs::NuisanceCorrelationsIndicator, m::EFTfitterModel, weights)

For all other models (i.e., covariance matrix is changed during sampling), this function 
computes the cross matrix without including observable weights in the matrix.
"""
function get_crossmatrix(mus::ModelUncertaintiesIndicator, ncs::NuisanceCorrelationsIndicator, m::EFTfitterModel, weights)
    cov = get_total_covariance(m)
    M_cov = m.CovarianceType(cov)

    return CovMatrix(M_cov, diag(cov))
end


# EFTfitterDensity constructor
function EFTfitterDensity(m::EFTfitterModel)
    measured_values = Float64[meas.value for meas in m.measurements]

    observable_functions = Function[meas.observable.prediction for meas in m.measurements]
    observable_mins = Float64[meas.observable.min for meas in m.measurements]
    observable_maxs = Float64[meas.observable.max for meas in m.measurements]
    weights = Float64[meas.observable.weight for meas in m.measurements]
    observable_weights = length(weights) * normalize(weights, 1)

    upper_bounds = any(x->x!=Inf, observable_maxs)
    lower_bounds = any(x->x!=-Inf, observable_mins)
    check_bounds = any([upper_bounds, lower_bounds]) ? BoundsCheck() : NoBoundsCheck()

    # check if model uncertainties are present
    #TODO: make this a function
    #TODO: add a warning if functions are too slow or take too much memory
    v = rand(m.parameters)
    predicted_values = [f(v) for f in observable_functions]
    mus = get_model_uncertainties_indicator(predicted_values)

    # preallocate arrays for storing parameter-dependent predictions and uncertainties 
    nthreads = Threads.nthreads()
    predictions = zeros(nthreads, length(observable_functions))
    prediction_uncertainties = zeros(nthreads, length(observable_functions))

    nuisance_correlations = build_nuisance_correlations(m.nuisances, m)

    #TODO: make this a function returning a matrix. depending on if there are model uncertainties or not, and if we want to use cholesky 
    crossmatrix = get_crossmatrix(mus, nuisance_correlations, m, weights)

    #TODO: add  support for limits
    limits = NoLimits()
    # limit_functions =  Function[]
    # limit_distributions = Distribution[]
    # limit_predictions = zeros(nthreads, length(limit_functions))
    # limit_uncertainties = zeros(nthreads, length(limit_functions))

    return EFTfitterDensity(
            measured_values,
            observable_functions,
            observable_mins,
            observable_maxs,
            observable_weights,
            crossmatrix,
            check_bounds,
            predictions,
            prediction_uncertainties,
            # limit_distributions,
            # limit_functions,
            # limit_predictions,
            # limit_uncertainties,
            mus,
            limits,
            nuisance_correlations
        )
end


function build_nuisance_correlations(nucs::Nothing, m::EFTfitterModel)
    return NoNuissanceCorrelations()
end 

function build_nuisance_correlations(nucs, m::EFTfitterModel)
    meas_keys = collect(keys(m.measurements))
    unc_keys = collect(keys(m.correlations))

    nuisances = NuisanceCorrelationIndices[]

    for (nui, nui_k) in zip(m.nuisances, collect(keys(m.nuisances)))
        println("nui: ", nui)
        println("nui_k: ", nui_k)
        unc = findfirst(x->x==nui.unc_key , unc_keys)
        i = findfirst(x->x==nui.meas1 , meas_keys)
        j = findfirst(x->x==nui.meas2 , meas_keys)
        push!(nuisances, NuisanceCorrelationIndices(unc, i, j, nui_k))
    end

    return  NuisanceCorrelations(nuisances, m)
end 



# build BAT.PosteriorMeasure directly from EFTfitterModel
function BAT.PosteriorMeasure(m::EFTfitterModel)
    likelihood = EFTfitterDensity(m)
    return posterior = BAT.PosteriorMeasure(likelihood, m.parameters)
end
