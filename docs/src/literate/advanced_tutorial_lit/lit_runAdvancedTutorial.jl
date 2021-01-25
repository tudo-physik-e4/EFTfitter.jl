#!jl # # EFTfitter.jl - Advanced Tutorial
#jl # EFTfitter.jl - Advanced Tutorial

#jl #~This tutorial introduces some of the advanced functionalities of EFTfitter.jl 
#jl #~using a generic example. Please see the [tutorial] for basic usage of EFTfitter.  
#md # Table of contents:
#md # ```@contents
#md # Pages = ["advanced_tutorial.md"]
#md # Depth = 3
#md # ```
                                                                                  #src
#nb using EFTfitter
#nb using BAT           # for sampling
#nb using IntervalSets  # for specifying the prior
#nb using Distributions # for specifying the prior
#nb using Plots         # for plotting
                                                                                 #src
#!jl include("lit_advanced_tutorial_inputs.jl")
#nb # ## File "runAdvancedTutorial.jl"
                                                                                 #src
#jl using EFTfitter
#jl using BAT            # for sampling
#jl using IntervalSets   # for specifying the prior
#jl using Distributions  # for specifying the prior
#jl using Plots # for plotting
#jl  
#jl include("advanced_inputs.jl")

#~We need to modify the definition of the `EFTfitterModel` by also passing the `nuisance_correlations`:
model = EFTfitterModel(parameters, measurements, correlations, nuisance_correlations)

#!md posterior = PosteriorDensity(model)
#!md 
#!md algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)
#!md samples = bat_sample(posterior, algorithm).result;

#!md p = plot(samples)
#!nb savefig(p, "plot.pdf")
#md # ![plot with nuisances](plots/plot_nuisance.png)
#md  # !!! warning
#md  #     When using nuisance correlations, it may happen that the fit does not converge anymore. 
#md  #     This happens mainly when the correlation values are close to -1 or +1. 
#md  #     In such a case, you can try to reduce the allowed ranges in the priors for the Nuisance 
#md  #     correlations to investigate at which values this happens.

#!md #~Warning: 
#!md #~When using nuisance correlations, it may happen that the fit does not converge anymore. 
#!md #~This happens mainly when the correlation values are close to -1 or +1. 
#!md #~In such a case, you can try to reduce the allowed ranges in the priors for the Nuisance 
#!md #~correlations to investigate at which values this happens.


#nb # ### Ranking of measurements and uncertainties
#md # ## Ranking of measurements and uncertainties
#jl #~============= Ranking of measurements and uncertainties ========================# 
#~With the `rank_measurements` and `rank_uncertainties` functions, the influence of
#~the individual measurements or uncertainty types on the result of a fit can be estimated.
#~For the ranking, each active measurement (respectively uncertainty type) is deactivated 
#~at a time and the fit is repeated. The results of all fits with one deactivated 
#~measurement (or uncertainty type) are then compared to the fit result with all measurements (uncertainty types) activated.
#~A ranking is calculated based on a ranking criterion calculated from the posterior distributions of these fits.

#~The default ranking criterion is the relative increase of the total width of the 
#~smallest interval containing 90% of the posterior probability when deactivating a measurement.
#~For models with more than one parameter, the sum of relative increases of all 
#~one-dimensional smallest intervals is used. The default ranking criterion is therefore
#`SumOfSmallestIntervals(p=0.9, bins=200)`.
measurement_ranking = EFTfitter.rank_measurements(model)

#~For ranking the uncertainty types, the relative decrease is used.
uncertainty_ranking = EFTfitter.rank_uncertainties(model, criterion = SumOfSmallestIntervals(p=0.9, bins=200))

#~Please see the [ranking documentation] for further ranking criteria and keyword arguments.
