export EFTfitterDensity
export EFTfitterDensityNuisance

struct _NuisanceCorrelation
    unc::Int
    i::Int
    j::Int
    key::Symbol # parameter key
end

struct EFTfitterDensity <: AbstractDensity
    measured_values::Vector{Float64}
    observable_functions::Vector{Function}
    observable_mins::Vector{Float64}
    observable_maxs::Vector{Float64}
    invcov::Array{Float64, 2}
    check_bounds::Bool
end

struct EFTfitterDensityNuisance <: AbstractDensity
    measured_values::Vector{Float64}
    observable_functions::Vector{Function}
    observable_mins::Vector{Float64}
    observable_maxs::Vector{Float64}
    covs::Vector{Array{Float64, 2}}
    nuisances::Vector{_NuisanceCorrelation}
    check_bounds::Bool
end


function EFTfitterDensity(m::EFTfitterModel)
    n = length(m.measurements)
    measured_values = [meas.value for meas in m.measurements]
    observable_functions = [meas.observable.func for meas in m.measurements]
    observable_mins = [meas.observable.min for meas in m.measurements]
    observable_maxs = [meas.observable.max for meas in m.measurements]

    bu = any(x->x!=Inf, observable_maxs)
    bl = any(x->x!=-Inf, observable_mins)
    check_bounds = any([bu, bl])


    invcov = inv(get_total_covariance(m))

    return EFTfitterDensity(
            measured_values,
            observable_functions,
            observable_mins,
            observable_maxs,
            invcov,
            check_bounds
            )
end

function EFTfitterDensityNuisance(m::EFTfitterModel)
    n = length(m.measurements)
    measured_values = [meas.value for meas in m.measurements]
    observable_functions = [meas.observable.func for meas in m.measurements]
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
    return all(withinbounds)
end


function evaluate_funcs(arr::Vector{Function}, params)
    return [arr[i](params) for i in 1:length(arr)]
end


function BAT.eval_logval_unchecked(
    m::EFTfitterDensity,
    params
)
    r = evaluate_funcs(m.observable_functions, params)

    if m.check_bounds
        ib = check_obs_bounds(r, m.observable_mins, m.observable_maxs)
        if ib == false
            return -Inf
        end
    end

    r = r-m.measured_values
    r1 = m.invcov*r
    result = -dot(r, r1)

    return  result
end


function BAT.eval_logval_unchecked(
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

    return  result
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



function BAT.PosteriorDensity(m::EFTfitterModel)
    if has_nuisance_correlations(m)
        likelihood = EFTfitterDensityNuisance(m)
        return posterior = PosteriorDensity(likelihood, m.parameters)
    else
        likelihood = EFTfitterDensity(m)
        return posterior = PosteriorDensity(likelihood, m.parameters)
    end
end
