export EFTfitterDensity
export EFTfitterDensityNuisance

struct _NuisanceCorrelation
    unc::Int
    i::Int
    j::Int
    key::Symbol # parameter key
end




struct EFTfitterDensity{M}
    measured_values::Vector{Float64}
    observable_functions::Vector{Function}
    observable_mins::Vector{Float64}
    observable_maxs::Vector{Float64}
    invcov::M
    check_bounds::Bool
    predictions::Matrix{Float64}
    prediction_uncertainties::Matrix{Float64}
    limit_distributions::Vector{Distribution}
    limit_functions::Vector{Function}
    limit_predictions::Matrix{Float64}
    limit_uncertainties::Matrix{Float64}
    has_model_uncertainties::Bool
    has_limits::Bool
    has_nuisance_correlations::Bool
end
@inline DensityInterface.DensityKind(::EFTfitterDensity) = IsDensity()



struct EFTfitterDensityNuisance
    measured_values::Vector{Float64}
    observable_functions::Vector{Function}
    observable_mins::Vector{Float64}
    observable_maxs::Vector{Float64}
    covs::Vector{Array{Float64, 2}}
    nuisances::Vector{_NuisanceCorrelation}
    check_bounds::Bool
end
@inline DensityInterface.DensityKind(::EFTfitterDensityNuisance) = IsDensity()


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


function EFTfitterDensity(m::EFTfitterModel)
    measured_values = Float64[meas.value for meas in m.measurements]
    observable_functions = Function[meas.observable.prediction for meas in m.measurements]
    observable_mins = Float64[meas.observable.min for meas in m.measurements]
    observable_maxs = Float64[meas.observable.max for meas in m.measurements]
    observable_weights = Float64[meas.observable.weight for meas in m.measurements]

    upper_bounds = any(x->x!=Inf, observable_maxs)
    lower_bounds = any(x->x!=-Inf, observable_mins)
    check_bounds = any([upper_bounds, lower_bounds])

    #todo: make this a function
    v = rand(m.parameters)
    predicted_values = [f(v) for f in observable_functions]
    has_model_uncertainties = typeof(predicted_values) == Vector{Float64} ? false : true

    invcov = inv(get_total_covariance(m))
    w = length(observable_weights) * normalize(observable_weights, 1)
    invcov_weighted = w .* invcov
    M_invcov = m.CovarianceType(invcov_weighted)

    nthreads = Threads.nthreads()
    predictions = zeros(nthreads, length(observable_functions))
    prediction_uncertainties = zeros(nthreads, length(observable_functions))

    has_limits = false
    limit_functions =  Function[]
    limit_distributions = Distribution[]
    limit_predictions = zeros(nthreads, length(limit_functions))
    limit_uncertainties = zeros(nthreads, length(limit_functions))

    has_nuisance_correlations = false

    return EFTfitterDensity(
            measured_values,
            observable_functions,
            observable_mins,
            observable_maxs,
            M_invcov,
            check_bounds,
            predictions,
            prediction_uncertainties,
            limit_distributions,
            limit_functions,
            limit_predictions,
            limit_uncertainties,
            has_model_uncertainties,
            has_limits,
            has_nuisance_correlations
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

function EFTfitterDensityNuisance(m::EFTfitterModel)
    n = length(m.measurements)
    measured_values = [meas.value for meas in m.measurements]
    observable_functions = [meas.observable.prediction for meas in m.measurements]
    observable_mins = [meas.observable.min for meas in m.measurements]
    observable_maxs = [meas.observable.max for meas in m.measurements]

    bu = any(x->x!=Inf, observable_maxs)
    bl = any(x->x!=-Inf, observable_mins)
    check_bounds = any([bu, bl])


    covs = get_covariances(m)

    meas_keys = collect(keys(m.measurements))
    unc_keys = collect(keys(m.correlations))


    nuisances = _NuisanceCorrelation[]
    for (nui, nui_k) in zip(m.nuisances, collect(keys(m.nuisances)))
        unc = findfirst(x->x==nui.unc_key , unc_keys)
        i = findfirst(x->x==nui.meas1 , meas_keys)
        j = findfirst(x->x==nui.meas2 , meas_keys)
        push!(nuisances, _NuisanceCorrelation(unc, i, j, nui_k))
    end

    return EFTfitterDensityNuisance(
            measured_values,
            observable_functions,
            observable_mins,
            observable_maxs,
            covs,
            nuisances,
            check_bounds
            )
end



function iswithinbounds(r::Float64, min::Float64, max::Float64)
    return min <= r <= max
end


function check_obs_bounds(r::Vector{Float64}, mins::Vector{Float64}, maxs::Vector{Float64})
    withinbounds = [iswithinbounds(r[i], mins[i], maxs[i]) for i in 1:length(r)]
    all(withinbounds) ? (return 1.) : (return -Inf)
end

# old:
# function evaluate_funcs(mu::Val{false}, arr::Vector{Function}, params)
#     return [arr[i](params) for i in 1:length(arr)]
# end

# without model uncertainties
function evaluate_funcs!(mu::Val{false}, arr::Vector{Function}, params, m)
    for i in eachindex(arr)
        m.predictions[Threads.threadid(), i] = arr[i](params)
    end
end

# with model uncertainties
function evaluate_funcs!(mu::Val{true}, arr::Vector{Function}, params, m)
    for i in eachindex(arr)
        res::Prediction = Prediction(arr[i](params))
        m.predictions[Threads.threadid(), i] = res.pred
        m.prediction_uncertainties[Threads.threadid(), i] = res.unc
    end
end

function add_uncs_to_invcov!(M, r)
    for i in 1:size(M)[1]
        M[i, i] += r[i]^2
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




# no nuisances and no model uncertainties
function DensityInterface.logdensityof(
    m::EFTfitterDensity,
    params
)
    evaluate_funcs!(Val(m.has_model_uncertainties), m.observable_functions, params, m)

    # #TODO: make a function
    # if m.check_bounds
    #     ib = check_obs_bounds(predictions, m.observable_mins, m.observable_maxs)
    #     if ib == false
    #         return -Inf
    #     end
    # end

    # TODO: change from dispatch on Vals to Types for speed
    result = calculate_l(m, Val(m.has_model_uncertainties), Val(m.has_limits))

    return result
end

# no model uncertainties, no limits
function calculate_l(m::EFTfitterDensity, mu::Val{false}, lim::Val{false})
    predictions = m.predictions[Threads.threadid(), :]

    inbounds = check_obs_bounds(predictions, m.observable_mins, m.observable_maxs)

    r = predictions-m.measured_values
    r1 = m.invcov*r
    result = -dot(r, r1)

    return  0.5*result*inbounds
end

# with model uncertainties
function calculate_l(m::EFTfitterDensity, mu::Val{true}, lim::Val{false})
    predictions = m.predictions[Threads.threadid(), :]

    inbounds = check_obs_bounds(predictions, m.observable_mins, m.observable_maxs)

    r = predictions-m.measured_values
    r1 = m.invcov*r
    result = -dot(r, r1)

    return  0.5*result*inbounds
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


function DensityInterface.logdensityof(
    m::EFTfitterDensityNuisance,
    params
)
    r = evaluate_funcs(m.observable_functions, params)

    if m.check_bounds
        ib = check_obs_bounds(r, m.observable_mins, m.observable_maxs)
        if ib == false
            return -Inf
        end
    end

    invcov = get_current_invcov(m, params)

    r = r-m.measured_values
    r1 = invcov*r
    result = -dot(r, r1)

    return  0.5*result
end


function get_current_invcov(m, params)
    for nui in m.nuisances
        i = nui.i; j = nui.j

        cov = params[nui.key] * sqrt(m.covs[nui.unc][i, i]) * sqrt(m.covs[nui.unc][j, j])
        m.covs[nui.unc][i, j] = cov
        m.covs[nui.unc][j, i] = cov
    end

    total_cov = sum(m.covs)
    invcov = inv(total_cov)
end



function BAT.PosteriorMeasure(m::EFTfitterModel)
    if has_nuisance_correlations(m)
        likelihood = EFTfitterDensityNuisance(m)
        return posterior = BAT.PosteriorMeasure(likelihood, m.parameters)
    else
        likelihood = EFTfitterDensity(m)
        return posterior = BAT.PosteriorMeasure(likelihood, m.parameters)
    end
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





    