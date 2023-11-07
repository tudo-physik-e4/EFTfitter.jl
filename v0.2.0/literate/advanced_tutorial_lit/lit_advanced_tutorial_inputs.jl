#jl # EFTfitter.jl - Advanced Tutorial

#jl #~This tutorial introduces some of the advanced functionalities of EFTfitter.jl 
#jl #~using a generic example. Please see the [tutorial] for basic usage of EFTfitter.
#nb # ---
#nb # It is recommended to store the inputs (i.e. parameters, observables, measurements and correlations) 
#nb # and the actual analysis in two separate files. This allows to simply load different configurations 
#nb # of the model. We will therefore also consider two files: `advanced_tutorial_inputs.jl` and `runAdvancedTutorial.jl`.

#nb # ## File "advanced_tutorial_inputs.jl"  

#nb # ### Parameters 
#jl #~============= Parameters =============================================# 
#!md #~We use the same parameters & priors as in the basic tutorial: 
#!md parameters = BAT.distprod(
#!md     C1 = -3..3, # short for: Uniform(-3, 3)
#!md     C2 = Normal(0, 0.5) # Normal distribution
#!md )


#nb # ### Observables
#jl #~============= Observables =============================================#
#!md #~We use the same observables as in the basic tutorial. However, we use a different
#!md #~way for creating the vector of functions for the BinnedMeasurement.

#!md function xsec1(params)
#!md     coeffs = [20.12, 5.56, 325.556]
#!md     return myfunc(params, coeffs)
#!md end
#!md 
#!md function xsec2(params)
#!md     coeffs = [2.12, 4.3, 12.6]
#!md     return myfunc(params, coeffs)
#!md end
#!md 
#!md function myfunc(params, c)
#!md     return c[1] * params.C1 + c[2] * params.C1 * params.C2+ c[3] * params.C2
#!md end


#md # ## Model uncertainties
#~The predictions for the observables can be affected by uncertainties that 
#~dependend on the current value of the model parameters. Such uncertainties can be included 
#~in the analysis by letting the observable functions return a tuple of numbers, where
#~the first number is the prediction and the second number is the absolute uncertainty on that prediction.
function xsec1(params)
    coeffs = [20.12, 5.56, 325.556]
    prediction = myfunc(params, coeffs)
    uncertainty = 0.1 * prediction # assume a 10% uncertainty
    return (prediction, uncertainty)
end

#~Note: It is possible to specify a model uncertainty only for some of the observables.
#~When the observable functions returns a single number, EFTfitter.jl will assume that the
#~model uncertainty is zero for this observable. However, please note that as soon as one observable function 
#~returns a tuple, the runtime of the model will probably increase, as the model uncertainties then need to be evaluate in each step.
#~Therefore, model uncertainties should only be used when the uncertainty value is depending on the
#~model parameters. Constant uncertainties should just be specified in the `uncertainties` field of the measurements. 



#-
#md # ## Vector of functions for a BinnedMeasurement
#~When using binned measurements, a vector of functions giving the predictions 
#~for the observable needs to be passed. It contains a function for each of bin and 
#~has only the model parameters as its argument. Defining a separate function for each 
#~bin can, however, become tedious for a large number of bins, especially since typically 
#~the bins of a distribution have a similar functional dependence on the model parameters 
#~and only differ in some coefficients. In such cases, it is possible to use Julia's 
#~anonymous functions to quickly create the vector of functions.  
#~The distribution in our [basic tutorial](https://tudo-physik-e4.github.io/EFTfitter.jl/dev/tutorial/) 
#~has been defined by implementing three functions that all call the same function `myfunc` 
#~but with different values for the coefficients
#~The same result can also be achieved like this:  

function get_coeffs(i) # return the coefficients for bin i
    coeffs = [[2.2, 5.5, 6.6], [2.2, 5.5, 6.6], [2.2, 5.5, 6.6]]
    return coeffs[i]
end 

function my_dist_func(params, i)
    coeffs = get_coeffs(i)
    return coeffs[1] * params.C1 + coeffs[2] * params.C1 * params.C2+ coeffs[3] * params.C2
end


# create an array of anonymous functions
diff_xsec = Function[x -> my_dist_func(x, i) for i in 1:3]

#nb # ### Measurements
#jl #~============= Measurements =============================================# 
#md # ## Using covariance matrices
#~Information about the uncertainties of measurements need to be provided to EFTfitter.jl 
#~in terms of the uncertainty values and corresponding correlation matrices.
#~If you have these information in terms of covariance matrices, you need to convert it 
#~to correlation matrices and uncertainty values before. The function 
#~[`cov_to_cor`](https://tudo-physik-e4.github.io/EFTfitter.jl/dev/api/#EFTfitter.cov_to_cor-Tuple{Array{var%22#s58%22,2}%20where%20var%22#s58%22%3C:Real}) 
#~can be used for this:

cov_syst = [3.24   0.81   0.378  0.324  0.468;
            0.81   0.81   0.126  0.162  0.234;
            0.378  0.126  0.49   0.126  0.182;
            0.324  0.162  0.126  0.81   0.234;
            0.468  0.234  0.182  0.234  1.69]

cor_syst, unc_syst = cov_to_cor(cov_syst)


#-
measurements = (

    Meas1 = Measurement(xsec1, 21.6,
            uncertainties = (stat=0.8, syst=unc_syst[1], another_unc=2.3)), 

    Meas2 = Measurement(Observable(xsec2, min=0), 1.9, 
            uncertainties = (stat=0.6, syst=unc_syst[2], another_unc=1.1), active=true),


    MeasDist = BinnedMeasurement(diff_xsec, [1.9, 2.93, 4.4],
                uncertainties = (stat = [0.7, 1.1, 1.2], syst= unc_syst[3:5], another_unc = [1.0, 1.2, 1.9]),
                active=[true, false, true]), 
)

#nb # ### Correlations
#jl #~============= Correlations =============================================# 
#!md #~The correlations are again the same as in the basic tutorial.
#!md dist_corr = [1.0 0.5 0.0;
#!md              0.5 1.0 0.0;
#!md              0.0 0.0 1.0]
#!md 
#!md another_corr_matrix = to_correlation_matrix(measurements,
#!md     (:Meas1, :Meas2, 0.4), 
#!md     (:Meas1, :MeasDist, 0.1), 
#!md     (:MeasDist, :MeasDist, dist_corr), 
#!md     (:MeasDist_bin2, :MeasDist_bin3, 0.3),
#!md )
#!md 
#!md #-
#!md 
#!md correlations = (
#!md     stat = NoCorrelation(active=true), # will use the identity matrix of the correct size
#!md 
#!md     syst = Correlation([1.0 0.5 0.3 0.2 0.2;
#!md                         0.5 1.0 0.2 0.2 0.2;
#!md                         0.3 0.2 1.0 0.2 0.2;
#!md                         0.2 0.2 0.2 1.0 0.2;
#!md                         0.2 0.2 0.2 0.2 1.0], active=true), # `active = false`: ignore all uncertainty values and correlations for this type of uncertainty
#!md 
#!md     another_unc = Correlation(another_corr_matrix, active=true) 
#!md )


#md # ## Nuisance Correlations
#nb # ### Nuisance Correlations
#jl #~============= Nuisance Correlations =============================================#
#~When performing an analysis with unknown correlation coefficients, it is possible 
#~to treat them as nuisance parameters in the fit.
#~For this, we define a further `NamedTuple` consisting of `NuisanceCorrelation` objects:

nuisance_correlations = (
    ρ1  = NuisanceCorrelation(:syst, :Meas1, :Meas2, -1..1),
    ρ2  = NuisanceCorrelation(:syst, :MeasDist_bin1, :MeasDist_bin3, truncated(Normal(0.5, 0.1), -1, 1)),
)

#~In the `NuisanceCorrelation` object we specify the name of the uncertainty type, the 
#~names of the two measurements we want to correlate using the nuisance correlations and 
#~a prior for the nuisance parameter. Note that the nuisance parameters should only be 
#~varied in the interval (-1, 1) as they represent correlation coefficients.
#~For `ρ1` we choose a flat prior between -1 and 1. For `ρ2` we have some expectations and
#~formulate them using a Gaussian prior with μ=0.5 and σ=0.1. However, to ensure that `ρ2` 
#~is only varied in the allowed region of (-1, 1), we [truncate](https://juliastats.org/Distributions.jl/stable/truncate/#Distributions.truncated) 
#~the normal distribution accordingly.
