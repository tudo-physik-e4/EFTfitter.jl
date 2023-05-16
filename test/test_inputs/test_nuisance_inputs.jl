using EFTfitter
using IntervalSets
using BAT

@testset "Test inputs with NuisanceCorrelation" begin
    
    parameters1 = BAT.NamedTupleDist(
        p1 = -20..20, 
        p2 = -10..10,
    )

    function testfunc1(params)
        c = [20.12, 5.56, 325.556]
        return c[1] * params.p1^2 + c[2] * params.p1 * params.p2 + c[3] * params.p2^2
    end

    measurements1 = (
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

    corr_matrix = to_correlation_matrix(measurements1,
        (:meas1, :meas2, 0.12), # will be ignored later in EFTfitterModel
        (:meas1, :meas3, 0.13), 
        (:meas1, :meas4_0_5, 0.141), 
        (:meas1, :meas4_5_10, 0.142), # will be ignored later in EFTfitterModel
        (:meas1, :meas4_10_20, 0.143), # will be ignored later in EFTfitterModel
        (:meas4_0_5, :meas4_5_10, 0.412), # will be ignored later in EFTfitterModel
        (:meas4_0_5, :meas4_5_10, 0.413), 
        (:meas4_0_5, :meas4_0_5, 0.9) # will be ignored later in EFTfitterModel
    )
    
    @test corr_matrix ≈ [ 1.0    0.12  0.13  0.141  0.142  0.143;
                         0.12   1.0   0.0   0.0    0.0    0.0;
                         0.13   0.0   1.0   0.0    0.0    0.0;
                         0.141  0.0   0.0   1.0    0.413  0.0;
                         0.142  0.0   0.0   0.413  1.0    0.0;
                         0.143  0.0   0.0   0.0    0.0    1.0] rtol = 0.01

    correlations1 = (
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

    model1 = EFTfitterModel(parameters1, measurements1, correlations1, limits=nothing, nuisances=nuisance_correlations, CovarianceType=Matrix)

    @testset "Test EFTfitterModel" begin
        @test length(model1.measurements) == 4
        @test model1.measurements.meas1 == Measurement(Observable(testfunc1, min=-Inf, max=Inf), 111.1, (unc1=10.1, unc3=13.3), true)
        @test model1.measurements.meas4_10_20 == Measurement(Observable(testfunc1, min=-Inf, max=Inf), 30., (unc1=0.13, unc3=0.33), true)
        @test keys(model1.measurements) == (:meas1, :meas3, :meas4_0_5, :meas4_10_20)  
        
        @test model1.measured_distributions.meas4.observable == [Observable(testfunc1), Observable(testfunc1)]
        @test model1.measured_distributions.meas4.value == [10., 30.]
        @test model1.measured_distributions.meas4.uncertainties == (unc1=[0.11, 0.13], unc3=[0.31, 0.33])
        @test model1.measured_distributions.meas4.active == [true, true]
        @test model1.measured_distributions.meas4.bin_names == [Symbol("0_5"), Symbol("10_20")] 
        
        @test length(model1.parameters._internal_distributions) == 3
        @test model1.parameters._internal_distributions.ρ1 == Uniform(0, 0.5)
        
        @test model1.nuisances == (ρ1=NuisanceCorrelation(:unc1, :meas1, :meas3, Uniform(0, 0.5)),)
    end
    
    @testset "Test EFTfitterDensity" begin
        eftfitter_density = PosteriorMeasure(model1).likelihood.density._d
        @test isa(eftfitter_density, EFTfitterDensityNuisance)
        
        @test eftfitter_density.measured_values == [111.1, 333.3, 10, 30]
        @test eftfitter_density.observable_functions == [testfunc1, testfunc1, testfunc1, testfunc1]
        @test eftfitter_density.observable_mins == [-Inf, 0, -Inf, -Inf]
        @test eftfitter_density.observable_maxs == [Inf, 1000.0, Inf, Inf]
        @test eftfitter_density.check_bounds == true
        @test eftfitter_density.covs ≈ [[102.01    0.0   0.0     0.0;
                                       0.0   906.01  0.0     0.0;
                                       0.0     0.0   0.0121  0.0;
                                       0.0     0.0   0.0     0.0169], 
                                      [176.89       52.3887  0.581343  0.627627;
                                      52.3887    918.09    0.0       0.0;
                                       0.581343    0.0     0.0961    0.0;
                                       0.627627    0.0     0.0       0.1089]] rtol=0.001
                                       
        @test eftfitter_density.nuisances == [EFTfitter._NuisanceCorrelation(1, 1, 2, :ρ1)]
    
        @test EFTfitter.get_current_invcov(eftfitter_density, (p1=0.0, p2=0.0, ρ1=0)) ≈ [0.00368812   -0.000105924  -0.0198158    -0.0184004;
                                                                                         -0.000105924   0.000551258   0.000569114   0.000528464;
                                                                                         -0.0198158     0.000569114   9.34861       0.0988625;
                                                                                         -0.0184004     0.000528464   0.0988625     8.04093] rtol=0.01
                                         
     @test EFTfitter.get_current_invcov(eftfitter_density, (p1=0.0, p2=0.0, ρ1=0.2)) ≈ [0.00376476   -0.000233614  -0.0202275   -0.0187827
                                                                                         -0.000233614   0.000562712   0.00125517   0.00116552;
                                                                                         -0.0202275     0.00125517    9.35082      0.100917;
                                                                                         -0.0187827     0.00116552    0.100917     8.04283] rtol=0.01                                 
        # test evaluation of correct likelihood value at a few points:
        @test DensityInterface.logdensityof(eftfitter_density, (p1=0.0, p2=0.0, ρ1=0)) ≈ -4088.80215 rtol=0.01
        @test DensityInterface.logdensityof(eftfitter_density, (p1=10.0, p2=-30.0, ρ1=0)) ≈ -Inf
        @test DensityInterface.logdensityof(eftfitter_density, (p1=1.5, p2=-0.9, ρ1=0)) ≈ -699198.52 rtol=0.01
        @test DensityInterface.logdensityof(eftfitter_density, (p1=-1.2, p2=-0.8, ρ1=0)) ≈ -438653.692 rtol=0.01
        
        @test DensityInterface.logdensityof(eftfitter_density, (p1=0.0, p2=0.0, ρ1=0.1)) ≈ -4091.12445 rtol=0.01
        @test DensityInterface.logdensityof(eftfitter_density, (p1=10.0, p2=-30.0, ρ1=0.2)) ≈ -Inf
        @test DensityInterface.logdensityof(eftfitter_density, (p1=1.5, p2=-0.9, ρ1=0.4)) ≈ -699985.42 rtol=0.01
    end
    
end
