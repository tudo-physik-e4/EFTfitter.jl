@testset "Test cov_to_cor" begin
        
    cov1 = [2.74 1.15 0.86 1.31;
            1.15 1.67 0.82 1.32;
            0.86 0.82 2.12 1.05;
            1.31 1.32 1.05 2.93]
            
    cor1, unc1 = EFTfitter.cov_to_cor(cov1)
    
    @test cor1 ≈[1.0 0.5376062677719541 0.35682497116489126 0.46234078289360164; 
                0.5376062677719541 1.0 0.4358004620316949 0.5967355144471027; 
                0.35682497116489126 0.4358004620316949 1.0 0.4212962160310593; 
                0.46234078289360164 0.5967355144471027 0.4212962160310593 1.0]
                
    @test unc1 ≈ [1.6552945357246849, 1.2922847983320085, 1.4560219778561037, 1.711724276862369]
    
    
    # test posdef warning
    cov2 = [2.74 -1.15 0.86 1.31;
            -1.15 1.67 0.82 1.32;
            0.86 0.82 2.12 1.05;
            1.31 1.32 1.05 2.93]
    @test_logs (:warn, "The covariance matrix [2.74 -1.15 0.86 1.31; -1.15 1.67 0.82 1.32; 0.86 0.82 2.12 1.05; 1.31 1.32 1.05 2.93] is not positive definite!") cor2, unc2 = EFTfitter.cov_to_cor(cov2)

end


@testset "Test all_observables_equal" begin
        estimator(params) = params.τ
        estimator2(params) = 2*params.τ
    
        parameters = BAT.NamedTupleDist(τ = 8..14, )

        covariance = [2.74 1.15 0.86 1.31;
                      1.15 1.67 0.82 1.32;
                      0.86 0.82 2.12 1.05;
                      1.31 1.32 1.05 2.93]

        corr, unc = EFTfitter.cov_to_cor(covariance)

        # all observables the same
        measurements = (
            τ1 = Measurement(estimator,  9.5, uncertainties = (stat=unc[1],) ),
            τ2 = Measurement(estimator, 11.9, uncertainties = (stat=unc[2],) ),
            τ3 = Measurement(estimator, 11.1, uncertainties = (stat=unc[3],) ),
            τ4 = Measurement(estimator,  8.9, uncertainties = (stat=unc[4],) ),
        )
        
        # not all observables the same
        measurements2 = (
            τ1 = Measurement(estimator,  9.5, uncertainties = (stat=unc[1],) ),
            τ2 = Measurement(estimator, 11.9, uncertainties = (stat=unc[2],) ),
            τ3 = Measurement(estimator2, 11.1, uncertainties = (stat=unc[3],) ),
            τ4 = Measurement(estimator,  8.9, uncertainties = (stat=unc[4],) ),
        )
        
        # all observables the same
        measurements3 = (
            τ1 = Measurement(estimator,  9.5, uncertainties = (stat=unc[1],) ),
            τ2 = Measurement(estimator, 11.9, uncertainties = (stat=unc[2],) ),
            τ3 = Measurement(estimator2, 11.1, uncertainties = (stat=unc[3],), active=false),
            τ4 = Measurement(estimator,  8.9, uncertainties = (stat=unc[4],) ),
        )
        
        correlations = (stat = Correlation(corr), )
        correlations2 = (stat = NoCorrelation(), )
        correlations3 = (stat = NoCorrelation(), )
        
        m = EFTfitterModel(parameters, measurements, correlations)
        m2 = EFTfitterModel(parameters, measurements2, correlations2)
        m3 = EFTfitterModel(parameters, measurements3, correlations3)
        
        @test  EFTfitter.all_observables_equal(m) == true
        @test  EFTfitter.all_observables_equal(m2) == false
        @test  EFTfitter.all_observables_equal(m3) == true
end


@testset "Test BLUE" begin
    estimator(params) = params.τ

    parameters = BAT.NamedTupleDist(τ = 8..14, )

    covariance = [2.74 1.15 0.86 1.31;
                  1.15 1.67 0.82 1.32;
                  0.86 0.82 2.12 1.05;
                  1.31 1.32 1.05 2.93]

    corr, unc = EFTfitter.cov_to_cor(covariance)

    measurements = (
        τ1 = Measurement(estimator,  9.5, uncertainties = (stat=unc[1],) ),
        τ2 = Measurement(estimator, 11.9, uncertainties = (stat=unc[2],) ),
        τ3 = Measurement(estimator, 11.1, uncertainties = (stat=unc[3],) ),
        τ4 = Measurement(estimator,  8.9, uncertainties = (stat=unc[4],) ),
    )
    
    correlations = (stat = Correlation(corr), )
    
    m = EFTfitterModel(parameters, measurements, correlations)
    
    blue = EFTfitter.BLUE(m)
    @test blue.value ≈ 11.1598305 rtol=0.005
    @test blue.unc ≈ 1.28604 rtol=0.005
    @test blue.weights ≈ [0.1450748, 0.46957738, 0.34729705, 0.0380508] rtol=0.005
end
