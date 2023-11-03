@testset "Test ranking" begin
    parameters1 = BAT.NamedTupleDist(
        p1 = -0.2..0.2, 
        p2 = -0.1..0.1,
    )

    function testfunc1(params)
    c = [20.12, 5.56, 325.556]
    return c[1] * params.p1^2 + c[2] * params.p1 * params.p2 + c[3] * params.p2^2
    end

    measurements1 = (
    meas1 = Measurement(testfunc1, 0.2, 
    uncertainties = (unc1=10.1, unc2=12.2, unc3=13.3), active=true),

    meas2 = Measurement(Observable(testfunc1, min=0, max=1000), 0.3, 
    uncertainties = (unc1=20.1, unc2=20.2, unc3=23.3), active=false),

    meas3 = Measurement(Observable(testfunc1, min=0, max=1000), 0.3, 
    uncertainties = (unc1=30.1, unc2=30.2, unc3=30.3), active=true),

    meas4 = BinnedMeasurement(Function[testfunc1, testfunc1, testfunc1],
    [0.28, 0.27, 0.29], uncertainties = (unc1=[0.11, 0.12, 0.13], unc2=[0.21, 0.22, 0.23], unc3=[0.31, 0.32, 0.33]),
    active = [true, false, true])
    )

    corr_matrix = to_correlation_matrix(measurements1,
    (:meas1, :meas2, 0.12), # will be ignored later in EFTfitterModel
    (:meas1, :meas3, 0.13), 
    (:meas1, :meas4_bin1, 0.141), 
    (:meas1, :meas4_bin2, 0.142), # will be ignored later in EFTfitterModel
    (:meas1, :meas4_bin3, 0.143), # will be ignored later in EFTfitterModel
    (:meas4_bin1, :meas4_bin2, 0.412), # will be ignored later in EFTfitterModel
    (:meas4_bin1, :meas4_bin2, 0.413), 
    (:meas4_bin1, :meas4_bin1, 0.9) # will be ignored later in EFTfitterModel
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

    model1 = EFTfitterModel(parameters1, measurements1, correlations1)
    
    measurement_models, measurement_keys = EFTfitter.measurement_models(model1)
    unc_models, unc_keys = EFTfitter.uncertainty_models(model1)
    
    @testset "Test model creation" begin
        measurement_models, measurement_keys = EFTfitter.measurement_models(model1)
        @test length(measurement_models) == 4
        @test measurement_keys == [:meas1, :meas3, :meas4_bin1, :meas4_bin3]
        @test keys(measurement_models[1].measurements) == (:meas3, :meas4_bin1, :meas4_bin3)
        @test keys(measurement_models[2].measurements) == (:meas1, :meas4_bin1, :meas4_bin3)
        
        
        unc_models, unc_keys = EFTfitter.uncertainty_models(model1)
        @test length(unc_models) == 2
        @test unc_keys == [:unc1, :unc3]
        @test keys(unc_models[1].correlations) == (:unc3,)
        @test keys(unc_models[2].correlations) == (:unc1,)
    end
    
    @testset "Test criterion_value" begin
        @test EFTfitter.criterion_value(SmallestInterval(key=:p1, p=0.9), unc_models[1],) ≈ 0.284 rtol = 0.05
        @test EFTfitter.criterion_value(SmallestInterval(key=:p2, p=0.9), unc_models[1],) ≈ 0.071 rtol = 0.05
        
        @test EFTfitter.criterion_value(SumOfSmallestIntervals(p=0.9), unc_models[1],) ≈ [0.284, 0.071] rtol = 0.05
    end
    
    @testset "Test ranking function" begin
        meas_ranks = EFTfitter.rank_measurements(model1, sampling_algorithm=SobolSampler(nsamples = 10^5))  
        @test meas_ranks.names == [:meas4_bin1, :meas4_bin3, :meas1, :meas3] 
        @test meas_ranks.values ≈ [0.229, 0.166, 0.0, 0.0] rtol=0.1

        unc_ranks = EFTfitter.rank_uncertainties(model1, sampling_algorithm=SobolSampler(nsamples = 10^5), order=:names, rev=false)  
        @test unc_ranks.names == [:unc1, :unc3] 
        @test unc_ranks.values ≈ [0.0278, 0.306] rtol=0.1
    end
end
