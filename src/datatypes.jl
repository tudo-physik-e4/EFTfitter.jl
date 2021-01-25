export Observable

export Measurement
export MeasurementDistribution

export Correlation
export NoCorrelation
export NuisanceCorrelation


"""
    struct Observable

Fields:  
* `func::Function`: Function returning the predicted value of the observable as a function of the parameters
* `min::Float64`: Minimum boundary for values of the observable. Defaults to `-Inf`.  
* `max::Float64`: Maximum boundary for values of the observable. Defaults to `Inf`.  

Constructors:

```julia
Observable(
    func::Function;
    min::Float64 = -Inf
    max::Float64 = Inf
)
```
"""
struct Observable
    func::Function
    min::Float64
    max::Float64
end

function Observable(f::Function; min=-Inf, max=Inf)
    Observable(f, min, max)
end


#----- Measurement -----------------------------------------
abstract type AbstractMeasurement end

"""
    struct Measurement

Fields:  
* `observable::Observable`: Observable that is measured.  
* `value::Float64;`: Measured value.   
* `uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Real}}}`: Uncertainties of the measurement as NamedTuple.  
* `active::Bool`: Use or exclude measurement in fit. Defaults to `true`.   

Constructors:
```julia
Measurement(
    observable::Observable,
    value::Float64;
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Real}}},
    active::Bool = true
)
```

```julia
Measurement(
    observable::Function,
    value::Float64;
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Real}}},
    active::Bool = true
)
```
"""
struct Measurement <: AbstractMeasurement
    observable::Observable
    value::Float64
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
    struct MeasurementDistribution

Fields:  
    * `observable::Array{Observable, 1}`: Observables that are measured.
    * `value::Array{Float64, 1}`: Measured values.
    * `uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Array{Float64, 1}}}}`: Uncertainties of the measurement as NamedTuple.
    * `active::Array{Bool, 1}`: Use or exclude bins in fit. Defaults to `true` for all bins.
    * `bin_names::Array{Symbol, 1}`: Suffixes that will be appended to the name of the measurement distribution for the individual bins. Defaults to [_bin1, _bin2, ...].

Constructors:  
```julia
MeasurementDistribution(
    observable::Array{Observable, 1},
    vals::Array{<:Real, 1};
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Union{Vector{Float64}, Vector{Int64}}}}},
    active::Union{Bool, Array{Bool, 1}} = [true for i in 1:length(vals)],
    bin_names::Array{Symbol, 1} = [Symbol("bin\$i") for i in 1:length(vals)]
)
```

```julia
MeasurementDistribution(
    observable::Array{Function, 1},
    vals::Array{<:Real, 1};
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Union{Vector{Float64}, Vector{Int64}}}}},
    active::Union{Bool, Array{Bool, 1}} = [true for i in 1:length(vals)],
    bin_names::Array{Symbol, 1}
)
```
"""
struct MeasurementDistribution <: AbstractMeasurement
    observable::Array{Observable, 1}
    value::Array{Float64, 1}
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Array{Float64, 1}}}}
    active::Array{Bool, 1}
    bin_names::Array{Symbol, 1}
end

# constructor with default value active=true
function MeasurementDistribution(
    observable::Array{Observable, 1},
    vals::Array{<:Real, 1};
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Union{Vector{Float64}, Vector{Int64}}}}},
    active::Union{Bool, Array{Bool, 1}} = [true for i in 1:length(vals)],
    bin_names::Array{Symbol, 1} = [Symbol("bin$i") for i in 1:length(vals)]
)
    isa(active, Bool) ? active = fill(active, length(vals)) : nothing
    unc = namedtuple(keys(uncertainties), float.(values(uncertainties)))

    MeasurementDistribution(observable, vals, unc, active, bin_names)
end

# constructor with default value active=true
function MeasurementDistribution(
    observable::Array{Function, 1},
    vals::Array{<:Real, 1};
    uncertainties::NamedTuple{<:Any, <:Tuple{Vararg{Union{Vector{Float64}, Vector{Int64}}}}},
    active::Union{Bool, Array{Bool, 1}} = [true for i in 1:length(vals)],
    bin_names::Array{Symbol, 1} = [Symbol("bin$i") for i in 1:length(vals)]
)
    obs = Observable.(observable)
    isa(active, Bool) ? active = fill(active, length(vals)) : nothing
    unc = namedtuple(keys(uncertainties), float.(values(uncertainties)))

    MeasurementDistribution(obs, vals, unc, active, bin_names)
end


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
struct Correlation <: AbstractCorrelation
    matrix::Array{Float64, 2}
    active::Bool
end

function Correlation(matrix::Array{<:Real, 2}; active::Bool = true)
    Correlation(matrix, active)
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
