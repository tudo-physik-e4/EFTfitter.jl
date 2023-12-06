export Observable
export Prediction

export Measurement
export BinnedMeasurement

export Correlation
export NoCorrelation
export NuisanceCorrelation


"""
    Prediction(pred::Float64, unc::Float64)

Represents the prediction, i.e. the nominal value and the associated absolute uncertainty, of an observable for a given set of parameters.

# Fields

- `pred::Float64`: The prediction value.
- `unc::Float64`: The absolute uncertainty associated with the prediction.

# Constructors

- `Prediction(a::Float64)`: Constructs a `Prediction` with a given value `a` and an uncertainty of `0.0`.
- `Prediction(a::Tuple{Float64, Float64})`: Constructs a `Prediction` from a tuple where the first element is the prediction value and the second is the uncertainty.
- `Prediction(a::Prediction)`: Identity constructor, returns the input `Prediction` object.
"""
struct Prediction
    pred::Float64
    unc::Float64
end

Prediction(a::Float64) = Prediction(a, 0.)
Prediction(a::Tuple{Float64, Float64}) = Prediction(a[1], a[2])
Prediction(a::Prediction) = a


"""
    struct Observable

Fields:  
* `prediction::Function`: Function returning the predicted value of the observable as a function of the parameters.
* `min::Real`: Minimum boundary for values of the observable. Defaults to `-Inf`.  
* `max::Real`: Maximum boundary for values of the observable. Defaults to `Inf`.  
* `weight::Real`: Weight of this observable in the combination. Defaults to `1.0`.  

Constructors:

```julia
Observable(
    prediction::Function
    min::Real = -Inf
    max::Real = Inf
    weight::Real = 1.0
)
```
"""
@with_kw struct Observable
    prediction::Function
    min::Real = -Inf
    max::Real = Inf
    weight::Real = 1.0
end

function Observable(prediction::Function; min=-Inf, max=Inf, weight=1.0)
    Observable(prediction=prediction, min=min, max=max, weight=weight)
end


#----- Measurement -----------------------------------------
abstract type AbstractMeasurement end

"""
    struct Measurement

Fields:  
* `observable::Observable`: Observable that is measured.  
* `value::Real;`: Measured value.   
* `uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Real}}}`: Uncertainties of the measurement as NamedTuple.  
* `active::Bool`: Use or exclude measurement in fit. Defaults to `true`.   

Constructors:
```julia
Measurement(
    observable::Observable,
    value::Real;
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Real}}},
    active::Bool = true
)
```

```julia
Measurement(
    observable::Function,
    value::Real;
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Real}}},
    active::Bool = true
)
```
"""
struct Measurement <: AbstractMeasurement
    observable::Observable
    value::Real
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Real}}} 
    active::Bool
end

function Measurement(
    observable::Observable,
    value::Real;
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Real}}},
    active::Bool = true
)
    Measurement(observable, value, uncertainties, active)
end

# Short-hand constructor that converts a function to an Observable object
function Measurement(
    observable::Function,
    value::Real;
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Real}}},
    active::Bool = true
)
    obs = Observable(observable)
    Measurement(obs, value, uncertainties, active)
end


#----- Measurement Distribution-----------------------------------------
"""
    struct BinnedMeasurement

Fields:  
    * `observable::Array{Observable, 1}`: Observables that are measured.
    * `values::Array{Float64, 1}`: Measured values.
    * `uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Array{Float64, 1}}}}`: Uncertainties of the measurement as NamedTuple.
    * `active::Array{Bool, 1}`: Use or exclude bins in fit. Defaults to `true` for all bins.
    * `bin_names::Array{Symbol, 1}`: Suffixes that will be appended to the name of the measurement distribution for the individual bins. Defaults to [_bin1, _bin2, ...].

Constructors:  
```julia
BinnedMeasurement(
    observable::Array{Observable, 1},
    values::Array{<:Real, 1};
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Union{Vector{Float64}, Vector{Int64}}}}},
    active::Union{Bool, Array{Bool, 1}} = [true for i in 1:length(vals)],
    bin_names::Array{Symbol, 1} = [Symbol("bin\$i") for i in 1:length(vals)]
)
```

```julia
BinnedMeasurement(
    observable::Array{Function, 1},
    values::Array{<:Real, 1};
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Union{Vector{Float64}, Vector{Int64}}}}},
    active::Union{Bool, Array{Bool, 1}} = [true for i in 1:length(vals)],
    bin_names::Array{Symbol, 1}
)
```
"""
struct BinnedMeasurement <: AbstractMeasurement
    observable::Array{Observable, 1}
    value::Array{Float64, 1}
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Array{Float64, 1}}}}
    active::Array{Bool, 1}
    bin_names::Array{Symbol, 1}
end

# constructor with default value active=true and names "bin1", ...
function BinnedMeasurement(
    observables::Array{Observable, 1},
    value::Array{<:Real, 1};
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Vector{<:Real}}}},
    active::Union{Bool, Array{Bool, 1}} = [true for i in eachindex(value)],
    bin_names::Array{Symbol, 1} = [Symbol("bin$i") for i in eachindex(value)]
)
    isa(active, Bool) ? active = fill(active, length(value)) : nothing
    unc = namedtuple(keys(uncertainties), float.(NamedTupleTools.values(uncertainties)))

    BinnedMeasurement(observables, value, unc, active, bin_names)
end

# constructor converting Function to Observable with default value active=true and names "bin1", ...
function BinnedMeasurement(
    observables::Array{Function, 1},
    value::Array{<:Real, 1};
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Vector{<:Real}}}},
    active::Union{Bool, Array{Bool, 1}} = [true for i in eachindex(value)],
    bin_names::Array{Symbol, 1} = [Symbol("bin$i") for i in eachindex(value)]
)
    obs = Observable.(observables)
    isa(active, Bool) ? active = fill(active, length(value)) : nothing
    unc = namedtuple(keys(uncertainties), float.(NamedTupleTools.values(uncertainties)))

    BinnedMeasurement(obs, value, unc, active, bin_names)
end

#----- TODO: Limits -----------------------------------------
abstract type AbstractLimit end

struct ExponentialUpperLimit <: AbstractLimit
    observable::Observable
    limit::Float64
    cl::Float64
    active::Bool
end

function ExponentialUpperLimit(
    observable::Function,
    limit::Float64,
    cl::Float64;
    active::Bool = true
)
    ExponentialUpperLimit(Observable(observable), limit, cl, active)
end

export ExponentialUpperLimit


struct GaussianUpperLimit <: AbstractLimit
    observable::Observable
    best_fit::Float64
    limit::Float64
    cl::Float64
    active::Bool
end

function GaussianUpperLimit(
    observable::Function,
    best_fit::Float64,
    limit::Float64,
    cl::Float64;
    active::Bool = true
)
    GaussianUpperLimit(Observable(observable), best_fit, limit, cl, active)
end

export GaussianUpperLimit

#----- Correlation -----------------------------------------
abstract type AbstractCorrelation end

"""
    struct Correlation

Fields:  
* `matrix::Array{Float64, 2}`: Observables that are measured.  
* `active::Bool`: Use this uncertainty category in fit. Defaults to `true`.  

Constructors:
```julia
Correlation(matrix::Array{<:Real, 2}; active::Bool = true)
```
"""
struct Correlation{M} <: AbstractCorrelation
    matrix::M #Symmetric{Float64,Array{Float64,2}}
    active::Bool
end

function Correlation(matrix::Array{<:Real, 2}; active::Bool = true, MatrixType = Symmetric)
    new_matrix = MatrixType(matrix)
    Correlation(new_matrix, active)
end

@with_kw struct NoCorrelation <: AbstractCorrelation
    active::Bool = true
end


"""
    struct NuisanceCorrelation

Fields:  
* `unc_key::Symbol`: Name of uncertainty category.  
* `meas1::Symbol`: Name of first measurement.  
* `meas2::Symbol`: Name of second measurement.
* `prior::Distribution`: Prior distribution. Accepts the type `Distribution` and all other 
                types accepted by BAT.NamedTupleDist, e.g. `Interval` or `Real`.
  

Constructors:  
```julia
NuisanceCorrelation(unc_key::Symbol, meas1::Symbol, meas2::Symbol, prior::Any)
```
"""
struct NuisanceCorrelation
    unc_key::Symbol
    meas1::Symbol
    meas2::Symbol
    prior::Distribution
end

function NuisanceCorrelation(unc_key::Symbol, meas1::Symbol, meas2::Symbol, prior::Any)
    dist, shape = ValueShapes._ntd_dist_and_shape(prior)
    return NuisanceCorrelation(unc_key, meas1, meas2, dist)
end
