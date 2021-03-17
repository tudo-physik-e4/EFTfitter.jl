export BLUE, cov_to_cor

"""
    cov_to_cor(cov::Array{<:Real, 2})
    
Convert a covariance matrix `cov` to a correlation matrix and a vector of uncertainty values.

Returns a matrix and a vector.
Throws a warning when the covariance matrix is not positive definite.

Example:

```julia
cor, unc = cov_to_cor(cov)
```
"""
function cov_to_cor(cov::Array{<:Real, 2})
    if !isposdef(cov)
        @warn "The covariance matrix $cov is not positive definite!"
    end
    
    cor = StatsBase.cov2cor(cov, sqrt.(diag(cov)))
    unc = sqrt.(diag(cov))
    
    return cor, unc
end


"""
    BLUE(model::EFTfitterModel)
    
Calculates the best linear unbiased estimator (BLUE) for multiple 
measurements of the same observable, according to https://www.sciencedirect.com/science/article/pii/0168900288900186.

Note: Works only for an `EFTfitterModel` where all measurements have the same observable.
If this is not the case, an error is thrown.

Returns a `NamedTuple` with the fields:
* `:value`: BLUE value
* `:unc`: BLUE uncertainty
* `:weights`: Array with the weights for each measurement


Example:

```julia
blue = BLUE(model)
println(blue.value, blue.unc, blue.weights)
```
"""
function BLUE(m::EFTfitterModel)
    if !all_observables_equal(m)
        @error "The measurements in the EFTfitterModel have different observables.
        Cannot calculate the best linear unbiased estimator (BLUE)."
        return
    end
    
    meas = [meas.value for meas in values(m.measurements)]
    nmeas = length(meas)
    covariance = get_total_covariance(m)
    
    u = ones(nmeas)
    α = inv(covariance)*u/(u'*inv(covariance)*u)
    τ_BLUE = dot(α, meas)
    σ_BLUE = sqrt(dot(α, covariance*α))
    
    return (value = τ_BLUE, unc = σ_BLUE, weights = α)
end

function all_observables_equal(model::EFTfitterModel)
    observable_functions = [m.observable.func for m in values(model.measurements)]
    all(y->y==observable_functions[1], observable_functions)
end
