function iswithinbounds(r::Float64, min::Float64, max::Float64)
    return min <= r <= max
end


function check_observable_bounds(cb::NoBoundsCheck, predictions, mins::Vector{Float64}, maxs::Vector{Float64})
    return 1.0
end

function check_observable_bounds(cb::BoundsCheck, predictions, mins::Vector{Float64}, maxs::Vector{Float64})
    r = view(predictions, Threads.threadid(), :)
    withinbounds = [iswithinbounds(r[i], mins[i], maxs[i]) for i in eachindex(r)]
    all(withinbounds) ? (return 1.) : (return 1e50)
end


# without model uncertainties 
function evaluate_funcs!(nomuncs::NoModelUncertainties, D::EFTfitterDensity, current_params)
    funcs = D.observable_functions
    for i in eachindex(funcs)
        res::Prediction = Prediction(funcs[i](current_params))
        D.predictions[Threads.threadid(), i] = res.pred
    end
end

# with model uncertainties
function evaluate_funcs!(muncs::HasModelUncertainties, D::EFTfitterDensity, current_params)
    funcs = D.observable_functions
    for i in eachindex(funcs)
        res::Prediction = Prediction(funcs[i](current_params))
        D.predictions[Threads.threadid(), i] = res.pred
        D.prediction_uncertainties[Threads.threadid(), i] = res.unc
    end
end


function DensityInterface.logdensityof(
    D::EFTfitterDensity,
    params
)
    evaluate_funcs!(D.model_uncertainties, D, params)

    bounds_factor = check_observable_bounds(D.check_bounds, D.predictions, D.observable_mins, D.observable_maxs)

    result = calculate_likelihood(D.model_uncertainties, D.nuisance_correlations, D, params)

    return bounds_factor * result
end

# no model uncertainties & no limits, weights are already included in the inverse covariance matrix
function calculate_likelihood(nomuncs::NoModelUncertainties, ns::NoNuissanceCorrelations, D::EFTfitterDensity, params)
    predictions = view(D.predictions, Threads.threadid(), :)

    @assert isa(D.crossmatrix, InvCovMatrix) #TODO: remove

    r = predictions-D.measured_values
    r1 = D.crossmatrix.m * r
    result = -dot(r, r1)

    return  0.5*result
end


# with model uncertainties, with nuisance correlations, no limits
function calculate_likelihood(muncs, nc, D::EFTfitterDensity, params)
    predictions = view(D.predictions, Threads.threadid(), :)
    prediction_uncertainties = view(D.prediction_uncertainties, Threads.threadid(), :)

    covmatrix = get_current_cov(nc, D, params)   

    add_model_uncertainties_to_cov!(muncs, covmatrix, prediction_uncertainties) 

    invcov = inv(covmatrix.m)

    r = D.observable_weights .* (predictions - D.measured_values) 
    r1 = invcov*r

    result = -dot(r, r1)

    return  0.5*result
end


# TODO: specify for cholesky lower ? (For larger matrices >100 cholesky is faster)
function add_model_uncertainties_to_cov!(mus::HasModelUncertainties, covmatrix, pred_unc)
    for i in eachindex(axes(covmatrix.m, 1))
        covmatrix.m[i, i] = covmatrix.diagonal[i] + pred_unc[i]^2 #TODO: check that original_diag is not affected by nuisance correlations
    end
end

function add_model_uncertainties_to_cov!(mus::NoModelUncertainties, covmatrix, pred_unc)
    nothing
end


function get_current_cov(nc::NuisanceCorrelations, D, params)
    for nui in D.nuisance_correlations.nuisance_correlations #TODO: rethink naming 
        i = nui.i; j = nui.j

        cov = params[nui.key] * sqrt(D.nuisance_correlations.covs[nui.unc][i, i]) * sqrt(D.nuisance_correlations.covs[nui.unc][j, j])
        D.nuisance_correlations.covs[nui.unc][i, j] = cov
        D.nuisance_correlations.covs[nui.unc][j, i] = cov
    end
   
    crossmatrix = CovMatrix(sum(D.nuisance_correlations.covs), D.crossmatrix.diagonal)
    return crossmatrix 
end


function get_current_cov(nc::NoNuissanceCorrelations, D, params)
    return D.crossmatrix
end
