@testset "Test API datatypes" begin
    test_function(x) = x

    # Test type Observable:
    obs1 = @inferred Observable(test_function)
    @test obs1.func == test_function
    @test obs1.min == -Inf
    @test obs1.max == Inf

    obs2 = @inferred Observable(test_function, min=-5, max=5)
    @test obs2.func == test_function
    @test obs2.min == -5
    @test obs2.max == 5

    # Test type Measurement:
    meas1 = @inferred Measurement(test_function, 50.0, uncertainties=(unc1=0.5, unc2=0.3), active=true)
    meas2 = @inferred Measurement(obs1, 50.0, uncertainties=(unc1=0.5, unc2=0.3), active=true)
    @test meas1 == meas2
    @test keys(meas1.uncertainties) == (:unc1, :unc2)

    #Test type MeasurementDistribution:
    measdist1 = MeasurementDistribution([obs1, obs1], [50.0, 49.9], 
                uncertainties=(unc1=[0.5, 0.4], unc2=[0.3, 0.29]), active=[true, false])
    measdist2 = MeasurementDistribution([obs1, obs1], [50.0, 49.9], 
                uncertainties=(unc1=[0.5, 0.4], unc2=[0.3, 0.29]), active=true)
    @test measdist1.active == [true, false]
    @test measdist2.active == [true, true]

    #Test type Correlation:
    corr1 = Correlation([1 0.5; 0.5 1])
    corr2 = Correlation([1 0.5; 0.5 1], active=false)
    @test corr1.active == true
    @test corr2.active == false
    @test corr1.matrix ==  [1 0.5; 0.5 1]
    
    #Test type NuisanceCorrelation:
    nui1 = @inferred NuisanceCorrelation(:stat, :meas1, :meas2, Normal(0, 1))
    @test nui1.unc_key == :stat
    @test nui1.meas1 == :meas1
    @test nui1.meas2 == :meas2
    @test nui1.prior == Normal(0, 1)
    
    nui2 = NuisanceCorrelation(:stat, :meas1, :meas2, -1..1)
    @test nui2.prior == Uniform(-1, 1)
    
    nui3 = NuisanceCorrelation(:stat, :meas1, :meas2, 0)
    @test nui3.prior == BAT.ConstValueDist(0)
end
