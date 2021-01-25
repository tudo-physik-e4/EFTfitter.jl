__precompile__(true)

module EFTfitter

using BAT
using Distributions
using LinearAlgebra
using StatsBase
using NamedTupleTools
using Setfield
using Parameters
using RecipesBase
using Requires
using ValueShapes

include("datatypes.jl")
include("EFTfitterModel.jl")
include("EFTfitterDensity.jl")
include("ranking/ranking.jl")
include("plotting/plotting.jl")
include("utils.jl")

const _PLOTS_MODULE = Ref{Union{Module,Nothing}}(nothing)
_plots_module() = _PLOTS_MODULE[]

function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" _PLOTS_MODULE[] = Plots
end

end # module
