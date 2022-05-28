using Test
using EFTfitter
using BAT
using Distributions
using IntervalSets
using DensityInterface

Test.@testset "Package EFTfitter" begin
    include("test_datatypes.jl")
    include("test_utils.jl")
    include("test_inputs/test_inputs.jl")
    include("test_ranking.jl")
end
