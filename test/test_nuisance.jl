using EFTfitter
using IntervalSets
using Distributions
using BAT
using DensityInterface
using BenchmarkTools
using Test

parameters = BAT.NamedTupleDist(
    p1 = -20..20, 
    p2 = -10..10,
)

function testfunc1(params)
    c = [20.12, 5.56, 325.556]
    return c[1] * params.p1^2 + c[2] * params.p1 * params.p2 + c[3] * params.p2^2
end

measurements = (
    meas1 = Measurement(testfunc1, 111.1, 
    uncertainties = (unc1=10.1, unc2=12.2, unc3=13.3), active=true),
    
    meas2 = Measurement(Observable(testfunc1, min=0, max=1000), 222.2, 
    uncertainties = (unc1=20.1, unc2=20.2, unc3=23.3), active=false),
    
    meas3 = Measurement(Observable(testfunc1, min=0, max=1000), 333.3, 
    uncertainties = (unc1=30.1, unc2=30.2, unc3=30.3), active=true),
    
    meas4 = MeasurementDistribution(Function[testfunc1, testfunc1, testfunc1],
    [10, 20, 30], uncertainties = (unc1=[0.11, 0.12, 0.13], unc2=[0.21, 0.22, 0.23], unc3=[0.31, 0.32, 0.33]),
    active = [true, false, true], bin_names=[Symbol("0_5"), Symbol("5_10"), Symbol("10_20")])
)

corr_matrix = to_correlation_matrix(measurements,
    (:meas1, :meas2, 0.12), # will be ignored later in EFTfitterModel
    (:meas1, :meas3, 0.13), 
    (:meas1, :meas4_0_5, 0.141), 
    (:meas1, :meas4_5_10, 0.142), # will be ignored later in EFTfitterModel
    (:meas1, :meas4_10_20, 0.143), # will be ignored later in EFTfitterModel
    (:meas4_0_5, :meas4_5_10, 0.412), # will be ignored later in EFTfitterModel
    (:meas4_0_5, :meas4_5_10, 0.413), 
    (:meas4_0_5, :meas4_0_5, 0.9) # will be ignored later in EFTfitterModel
)

correlations = (
    unc1 = NoCorrelation(active=true),

    # wrong matrix size for number of measurements, will be ignored if active=false:
    unc2 = Correlation([1.0 0.5 0.7; 
                        0.5 1.0 0.6;
                        0.7 0.6 1.0], active=false),
    
    unc3 = Correlation(corr_matrix)
)

nuisance_correlations = (
    ρ1 = NuisanceCorrelation(:unc1, :meas1, :meas3, 0..0.5),
    ρ2 = NuisanceCorrelation(:unc1, :meas1, :meas2, truncated(Normal(0, 1), 0, 0.9)),
)


model = EFTfitterModel(parameters, measurements, correlations, nuisances=nuisance_correlations, )
posterior = PosteriorMeasure(model)

v = (p1 = 10.826122384321511, p2 = -8.32129957354641, ρ1 = 0.4, ρ2 = 0.8)
logp = logdensityof(posterior)

@test logp(v) ≈ 

@btime logp(v)

t = @benchmark logp(v)
@test t.allocs == 15
@test t.memory == 720
@test minimum(t.times) ≈ 226 atol=5
