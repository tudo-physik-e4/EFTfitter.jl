export EFTfitterDensity
export EFTfitterDensityNuisance

struct _NuisanceCorrelation
    unc::Int
    i::Int
    j::Int
    key::Symbol # parameter key
end

abstract type MatrixType end

struct CovarianceMatrix{T} <: MatrixType 
    m::T
end
struct InverseCovarianceMatrix{T} <: MatrixType 
    m::T
end

abstract type ModelUncertaintiesStatus end

struct HasModelUncertainties <: ModelUncertaintiesStatus end
struct NoModelUncertainties <: ModelUncertaintiesStatus end

get_model_uncertainties_status(predictions::Vector{<:Real}) = NoModelUncertainties()
get_model_uncertainties_status(predictions) = HasModelUncertainties()

abstract type AbstractNuisanceCorrelations end
struct NoNuissanceCorrelations <: AbstractNuisanceCorrelations end

struct NuisanceCorrelations{T} <: AbstractNuisanceCorrelations 
    nuisances::Vector{_NuisanceCorrelation}
    covs::Vector{T}
end

NuisanceCorrelations(nuisances, m::EFTfitterModel) = NuisanceCorrelations(nuisances, get_covariances(m))



abstract type LimitsStatus end
struct HasLimits <: LimitsStatus end
struct NoLimits <: LimitsStatus end


struct EFTfitterDensity{M<:MatrixType, MU<:ModelUncertaintiesStatus, NC<:AbstractNuisanceCorrelations, L<:LimitsStatus}
    measured_values::Vector{Float64}
    observable_functions::Vector{Function}
    observable_mins::Vector{Float64}
    observable_maxs::Vector{Float64}
    weights::Vector{Float64}
    matrix::M
    original_diag::Vector{Float64}
    check_bounds::Bool
    predictions::Matrix{Float64}
    prediction_uncertainties::Matrix{Float64}
    limit_distributions::Vector{Distribution}
    limit_functions::Vector{Function}
    limit_predictions::Matrix{Float64}
    limit_uncertainties::Matrix{Float64}
    mus::MU
    ls::L
    nuisance_correlations::NC
end
@inline DensityInterface.DensityKind(::EFTfitterDensity) = IsDensity()


# TODO: remove

# struct EFTfitterDensityNuisance
#     measured_values::Vector{Float64}
#     observable_functions::Vector{Function}
#     observable_mins::Vector{Float64}
#     observable_maxs::Vector{Float64}
#     covs::Vector{Array{Float64, 2}}
#     nuisances::Vector{_NuisanceCorrelation}
#     check_bounds::Bool
# end
# @inline DensityInterface.DensityKind(::EFTfitterDensityNuisance) = IsDensity()


# struct EFTfitterDensityWithLimits
#     measured_values::Vector{Float64}
#     observable_functions::Vector{Function}
#     observable_mins::Vector{Float64}
#     observable_maxs::Vector{Float64}
#     limit_distributions::Vector{Distribution}
#     limit_functions::Vector{Function}
#     invcov::Array{Float64, 2}
#     check_bounds::Bool
# end
# @inline DensityInterface.DensityKind(::EFTfitterDensityWithLimits) = IsDensity()

function prepare_cov_matrix(mus::NoModelUncertainties, ncs::NoNuissanceCorrelations, m::EFTfitterModel, weights)
    invcov = inv(get_total_covariance(m))
    
    invcov_weighted = weights .* invcov
    M_invcov = m.CovarianceType(invcov_weighted)

    return InverseCovarianceMatrix(M_invcov), diag(invcov_weighted)
end

function prepare_cov_matrix(mus::ModelUncertaintiesStatus, ncs::AbstractNuisanceCorrelations, m::EFTfitterModel, weights)
    cov = get_total_covariance(m)
    M_cov = m.CovarianceType(cov)

    return CovarianceMatrix(M_cov), diag(cov)
end


function EFTfitterDensity(m::EFTfitterModel)
    measured_values = Float64[meas.value for meas in m.measurements]
    observable_functions = Function[meas.observable.prediction for meas in m.measurements]
    observable_mins = Float64[meas.observable.min for meas in m.measurements]
    observable_maxs = Float64[meas.observable.max for meas in m.measurements]

    observable_weights = Float64[meas.observable.weight for meas in m.measurements]
    weights = length(observable_weights) * normalize(observable_weights, 1)

    upper_bounds = any(x->x!=Inf, observable_maxs)
    lower_bounds = any(x->x!=-Inf, observable_mins)
    check_bounds = any([upper_bounds, lower_bounds])

    #todo: make this a function
    #TODO: add a warning if functions are too slow or take too much memory
    v = rand(m.parameters)
    predicted_values = [f(v) for f in observable_functions]
    mus = get_model_uncertainties_status(predicted_values)

    nthreads = Threads.nthreads()
    predictions = zeros(nthreads, length(observable_functions))
    prediction_uncertainties = zeros(nthreads, length(observable_functions))

    #TODO: add limits
    ls = NoLimits()
    limit_functions =  Function[]
    limit_distributions = Distribution[]
    limit_predictions = zeros(nthreads, length(limit_functions))
    limit_uncertainties = zeros(nthreads, length(limit_functions))

    #TODO: add nuisance correlations
    meas_keys = collect(keys(m.measurements))
    unc_keys = collect(keys(m.correlations))

    nuisances = _NuisanceCorrelation[]
    for (nui, nui_k) in zip(m.nuisances, collect(keys(m.nuisances)))
        unc = findfirst(x->x==nui.unc_key , unc_keys)
        i = findfirst(x->x==nui.meas1 , meas_keys)
        j = findfirst(x->x==nui.meas2 , meas_keys)
        push!(nuisances, _NuisanceCorrelation(unc, i, j, nui_k))
    end

    nui =  isempty(nuisances) ? NoNuissanceCorrelations() : NuisanceCorrelations(nuisances, m)

    #TODO: make this a function returning a matrix. depending on if there are model uncertainties or not, and if we want to use cholesky 
    matrix, original_diagonal = prepare_cov_matrix(mus, ncs, m, weights)


    return EFTfitterDensity(
            measured_values,
            observable_functions,
            observable_mins,
            observable_maxs,
            weights,
            matrix,
            original_diagonal,
            check_bounds,
            predictions,
            prediction_uncertainties,
            limit_distributions,
            limit_functions,
            limit_predictions,
            limit_uncertainties,
            mus,
            ls,
            nui
        )
end


using SpecialFunctions
function Normal_from_limit(best_fit_value, limit, confidence_level)
    μ = best_fit_value
    p = confidence_level
    q = limit

    σ = (q - μ)/(sqrt(2)*erfinv(2*p-1))

    return Normal(μ, σ)
end 

function make_dist(limit::GaussianUpperLimit)
    return Normal_from_limit(limit.best_fit, limit.limit, limit.cl)
end

# function EFTfitterDensityWithLimits(m::EFTfitterModelWithLimits)
#     n = length(m.measurements)
#     measured_values = [meas.value for meas in m.measurements]
#     observable_functions = [meas.observable.prediction for meas in m.measurements]
#     observable_mins = [meas.observable.min for meas in m.measurements]
#     observable_maxs = [meas.observable.max for meas in m.measurements]

#     bu = any(x->x!=Inf, observable_maxs)
#     bl = any(x->x!=-Inf, observable_mins)
#     check_bounds = any([bu, bl])

#     limit_observable_functions = [limit.observable.prediction for limit in m.limits]
#     limit_distributions = [make_dist(limit) for limit in m.limits]

#     invcov = inv(get_total_covariance(m))

#     return EFTfitterDensityWithLimits(
#             measured_values,
#             observable_functions,
#             observable_mins,
#             observable_maxs,
#             limit_distributions,
#             limit_observable_functions,
#             invcov,
#             check_bounds
#             )
# end

# function EFTfitterDensityNuisance(m::EFTfitterModel)
#     n = length(m.measurements)
#     measured_values = [meas.value for meas in m.measurements]
#     observable_functions = [meas.observable.prediction for meas in m.measurements]
#     observable_mins = [meas.observable.min for meas in m.measurements]
#     observable_maxs = [meas.observable.max for meas in m.measurements]

#     bu = any(x->x!=Inf, observable_maxs)
#     bl = any(x->x!=-Inf, observable_mins)
#     check_bounds = any([bu, bl])


#     covs = get_covariances(m)

#     meas_keys = collect(keys(m.measurements))
#     unc_keys = collect(keys(m.correlations))


#     nuisances = _NuisanceCorrelation[]
#     for (nui, nui_k) in zip(m.nuisances, collect(keys(m.nuisances)))
#         unc = findfirst(x->x==nui.unc_key , unc_keys)
#         i = findfirst(x->x==nui.meas1 , meas_keys)
#         j = findfirst(x->x==nui.meas2 , meas_keys)
#         push!(nuisances, _NuisanceCorrelation(unc, i, j, nui_k))
#     end

#     return EFTfitterDensityNuisance(
#             measured_values,
#             observable_functions,
#             observable_mins,
#             observable_maxs,
#             covs,
#             nuisances,
#             check_bounds
#             )
# end



function iswithinbounds(r::Float64, min::Float64, max::Float64)
    return min <= r <= max
end

# TODO: add dispatch, return 1 by default
function check_obs_bounds(r::Vector{Float64}, mins::Vector{Float64}, maxs::Vector{Float64})
    withinbounds = [iswithinbounds(r[i], mins[i], maxs[i]) for i in 1:length(r)]
    all(withinbounds) ? (return 1.) : (return -Inf)
end


# why is this slower?
# function evaluate_funcs!(mus::NoModelUncertainties, arr::Vector{Function}, params, m)
#     for i in eachindex(arr)
#         m.predictions[Threads.threadid(), i] = arr[i](params)
#     end
# end


# without model uncertainties 
function evaluate_funcs!(mus::NoModelUncertainties, arr::Vector{Function}, params, m)
    for i in eachindex(arr)
        res::Prediction = Prediction(arr[i](params))
        m.predictions[Threads.threadid(), i] = res.pred
    end
end

# with model uncertainties
function evaluate_funcs!(mus::HasModelUncertainties, arr::Vector{Function}, params, m)
    for i in eachindex(arr)
        res::Prediction = Prediction(arr[i](params))
        m.predictions[Threads.threadid(), i] = res.pred
        m.prediction_uncertainties[Threads.threadid(), i] = res.unc
    end
end



# # old:
# function DensityInterface.logdensityof(
#     m::EFTfitterDensity,
#     params
# )
#     r = evaluate_funcs(m.observable_functions, params)

#     if m.check_bounds
#         ib = check_obs_bounds(r, m.observable_mins, m.observable_maxs)
#         if ib == false
#             return -Inf
#         end
#     end

#     r = r-m.measured_values
#     r1 = m.invcov*r
#     result = -dot(r, r1)

#     return  0.5*result
# end



function DensityInterface.logdensityof(
    m::EFTfitterDensity,
    params
)
    #Todo: move into calculate_likelihood?
    evaluate_funcs!(m.mus, m.observable_functions, params, m)

    #TODO: check_observable_bounds


    result = calculate_likelihood(m, m.mus, m.ls)

    return result
end

# no model uncertainties, no limits, weights are already in the covariance matrix
function calculate_likelihood(m::EFTfitterDensity, mus::NoModelUncertainties, ls::NoLimits)
    predictions = view(m.predictions, Threads.threadid(), :)

    @assert isa(m.matrix, InverseCovarianceMatrix)

    r = predictions-m.measured_values
    r1 = m.matrix.m*r
    result = -dot(r, r1)

    return  0.5*result
end


# TODO: specify for cholesky lower
function add_model_uncertainties!(m, r)
    for i in 1:size(m.matrix.m, 1)
        m.matrix.m[i, i] = m.original_diag[i] + r[i]^2
    end
end

# with model uncertainties, no limits
function calculate_likelihood(m::EFTfitterDensity, mus::HasModelUncertainties, ls::NoLimits)
    predictions = view(m.predictions, Threads.threadid(), :)
    prediction_uncertainties = view(m.prediction_uncertainties, Threads.threadid(), :)

    @assert isa(m.matrix, CovarianceMatrix)
    add_model_uncertainties!(m, prediction_uncertainties)
    invcov = inv(m.matrix.m)

    r = m.weights .* (predictions-m.measured_values) # we
    r1 = invcov*r

    result = -dot(r, r1)

    return  0.5*result
end

# function DensityInterface.logdensityof(
#     m::EFTfitterDensityWithLimits,
#     params
# )
#     r = evaluate_funcs(m.observable_functions, params)

#     if m.check_bounds
#         ib = check_obs_bounds(r, m.observable_mins, m.observable_maxs)
#         if ib == false 
#             return -Inf
#         end
#     end

#     r = r-m.measured_values
#     r1 = m.invcov*r
#     result = -dot(r, r1)
#     gaussian_result = 0.5*result

    
#     r_limits = evaluate_funcs(m.limit_functions, params)
#     ls=0.
#     for i in eachindex(m.limit_distributions)
#         ls += pdf(m.limit_distributions[i], r_limits[i])
#     end


#     return gaussian_result + ls
# end


# function DensityInterface.logdensityof(
#     m::EFTfitterDensityNuisance,
#     params
# )
#     r = evaluate_funcs(m.observable_functions, params)

#     if m.check_bounds
#         ib = check_obs_bounds(r, m.observable_mins, m.observable_maxs)
#         if ib == false
#             return -Inf
#         end
#     end

#     invcov = get_current_invcov(m, params)

#     r = r-m.measured_values
#     r1 = invcov*r
#     result = -dot(r, r1)

#     return  0.5*result
# end


function get_current_invcov(m, params)
    for nui in m.nuisance_correlations.nuisances
        i = nui.i; j = nui.j

        cov = params[nui.key] * sqrt(m.nuisance_correlations.covs[nui.unc][i, i]) * sqrt(m.nuisance_correlations.covs[nui.unc][j, j])
        m.nuisance_correlations.covs[nui.unc][i, j] = cov
        m.nuisance_correlations.covs[nui.unc][j, i] = cov
    end

    total_cov = sum(m.nuisance_correlations.covs)
    invcov = inv(total_cov)
end



function BAT.PosteriorMeasure(m::EFTfitterModel)
    # if has_nuisance_correlations(m) # TODO: remove
    #     likelihood = EFTfitterDensityNuisance(m)
    #     return posterior = BAT.PosteriorMeasure(likelihood, m.parameters)
    # else
    #     likelihood = EFTfitterDensity(m)
    #     return posterior = BAT.PosteriorMeasure(likelihood, m.parameters)
    # end
    likelihood = EFTfitterDensity(m)
    return posterior = BAT.PosteriorMeasure(likelihood, m.parameters)
end

# function BAT.PosteriorMeasure(m::EFTfitterModelWithLimits)
#     likelihood = EFTfitterDensityWithLimits(m)
#     return posterior = BAT.PosteriorMeasure(likelihood, m.parameters)

# end


# x = rand(5)
# M = rand(5,5)
# w = [4,4,4,4,5]
# w = length(w)*normalize(w, 1)

# mean(w)


# x = ones(5)

# r = dot(w.*x, M*x)

# w .* x
# x

# M2 = w .* M
# r2 = dot(x, M2*x)

# r_unw = dot(x, M*x)

# r_unw/r

# w = [1, 2, 2, 1, 3]
# w = normalize(w, 1)

# M2 =x .* M

# x2 = w.*x
# r = dot(x2, M*x)

# r2 = dot(x, M2*x)





    