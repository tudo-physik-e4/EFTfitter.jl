@testset "Test simple inputs" begin
    
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
    )

    correlations1 = (
        unc1 = NoCorrelation(active=true),

        unc2 = Correlation([1.0 0.5 0.7;
                            0.5 1.0 0.6;
                            0.7 0.6 1.0], active=false),
        
        unc3 = Correlation([1.0 0.2 0.3;
                            0.2 1.0 0.4;
                            0.3 0.4 1.0])
    )

    model1 = EFTfitterModel(parameters1, measurements1, correlations1)
    
    @testset "Test EFTfitterModel" begin
        @test model1.measurements.meas1.observable == Observable(testfunc1)

        @test length(model1.measurements) == 2
        @test keys(model1.measurements) == (:meas1, :meas3)
        @test model1.measurements.meas1 == Measurement(Observable(testfunc1, min=-Inf, max=Inf), 111.1, (unc1=10.1, unc3=13.3), true)
        @test model1.measurements.meas3 == Measurement(Observable(testfunc1, min=0, max=1000), 333.3, (unc1=30.1, unc3=30.3), true)        
    
        @test get_observables(model1) == (testfunc1 = Observable(testfunc1, min=-Inf, max=Inf),)
    end
    
    @testset "Test EFTfitterDensity" begin
        PosteriorMeasure(model1)
        eftfitter_density = PosteriorMeasure(model1).likelihood.density._d
        
        @test isa(eftfitter_density, EFTfitterDensity)
        
        @test eftfitter_density.measured_values == [111.1, 333.3]
        @test eftfitter_density.observable_functions == [testfunc1, testfunc1]
        @test eftfitter_density.observable_mins == [-Inf, 0]
        @test eftfitter_density.observable_maxs == [Inf, 1000.0]
        @test eftfitter_density.crossmatrix.m ≈ [ 0.00369157   -0.000244669; -0.000244669   0.000564432] rtol=0.01
        @test eftfitter_density.check_bounds == EFTfitter.BoundsCheck()
    
        # test evaluation of correct likelihood value at a few points:
        @test DensityInterface.logdensityof(eftfitter_density, (p1=0.0, p2=0.0)) ≈ -45.07398 rtol=0.01
        @test DensityInterface.logdensityof(eftfitter_density, (p1=10.0, p2=-30.0)) ≈ -1.6191906947178624e58  rtol=0.01
        @test DensityInterface.logdensityof(eftfitter_density, (p1=1.5, p2=-0.9)) ≈ -68.6575 rtol=0.01
        @test DensityInterface.logdensityof(eftfitter_density, (p1=-1.2, p2=-0.8)) ≈ -37.1857 rtol=0.01
    end
    
end

