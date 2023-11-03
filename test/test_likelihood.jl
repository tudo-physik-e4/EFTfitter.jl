using Random

using Random, LinearAlgebra

function generate_correlation_matrix(N)
    data = randn(N, N)  # Generate a random N-by-N matrix
    
    # Create a covariance matrix by multiplying the data matrix with its transpose
    covariance_matrix = data * data'
    
    # Calculate the diagonal matrix of standard deviations
    diagonal_matrix = diag(covariance_matrix)
    standard_deviations = sqrt.(diagonal_matrix)
    
    # Calculate the inverse of the diagonal matrix of standard deviations
    inv_diagonal_matrix = diagm(1 ./ standard_deviations)
    
    # Calculate the correlation matrix
    correlation_matrix = inv_diagonal_matrix * covariance_matrix * inv_diagonal_matrix
    
    return correlation_matrix
end



function create_model(N; use_model_uncertainties = false)
    Random.seed!(1234)

    parameters = BAT.NamedTupleDist(
        p1 = -20..20, 
        p2 = -10..10,
    )

    true_params = (p1=1.3, p2=2.5)

    function testfunc(params)
        c = [20.12, 5.56, 325.556]
        return c[1] * params.p1^2 + c[2] * params.p1 * params.p2 + c[3] * params.p2^2
    end

    function testfunc_mu(params)
        c = [20.12, 5.56, 325.556]
        result = c[1] * params.p1^2 + c[2] * params.p1 * params.p2 + c[3] * params.p2^2
        return (result, 0.1 )
    end


    expected = testfunc(true_params)

    
    func_arr = use_model_uncertainties ? Function[testfunc_mu for i in 1:N] : Function[testfunc for i in 1:N]
    
    meas_arr = [expected*abs(rand(Normal(1, 0.1))) for i in 1:N]
    unc_arr = [expected*abs(rand(Normal(0.1, 0.1))) for i in 1:N]


    measurements = (
        #meas1 = Measurement(testfunc, 3000.0, uncertainties = (unc = 50.1,)),

        meas = BinnedMeasurement(func_arr, meas_arr,
                    uncertainties = (unc = unc_arr,)),
    )

    corm = generate_correlation_matrix(N)

    correlations = (
        unc = Correlation(corm),
    )

    return EFTfitterModel(parameters, measurements, correlations)
end 

model = create_model(5, use_model_uncertainties = true)


@testset "Test likelihood" begin
    model_5 = create_model(5)
    posterior_5 = PosteriorMeasure(model_5)
    v = (p1 = 1.36, p2 = 2.3)
    @test logdensityof(posterior_5)(v) ≈ -39.5175576223029 rtol=0.005
   # @elapsed logdensityof(posterior_5)(v)


    model_5_mu = create_model(5, use_model_uncertainties = true)
    posterior_5_mu = PosteriorMeasure(model_5_mu)
    v = (p1 = 1.36, p2 = 2.3)
    @test logdensityof(posterior_5_mu)(v) ≈ -39.5175576223029 rtol=0.005



    model_500 = create_model(500)
    posterior_500 = PosteriorMeasure(model_500)
    v = (p1 = 1.36, p2 = 2.3)
    @test logdensityof(posterior_500)(v) ≈ -4.916679677894212e11 rtol=0.005
    #@btime logdensityof(posterior_500)(v)
    
    model_5000 = create_model(5000)
    posterior_5000 = PosteriorMeasure(model_5000)
    v = (p1 = 1.36, p2 = 2.3)
    @test logdensityof(posterior_5000)(v) ≈ -3.706920707573017e11 rtol=0.005
    #@btime logdensityof(posterior_5000)(v) 
end