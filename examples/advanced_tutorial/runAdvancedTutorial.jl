# EFTfitter.jl - Advanced Tutorial
# This tutorial introduces some of the advanced functionalities of EFTfitter.jl
# using a generic example. Please see the [tutorial] for basic usage of EFTfitter.
using EFTfitter
using BAT            # for sampling
using IntervalSets   # for specifying the prior
using Distributions  # for specifying the prior
using Plots # for plotting

include("advanced_inputs.jl")

# We need to modify the definition of the `EFTfitterModel` by also passing the `nuisance_correlations`:
model = EFTfitterModel(parameters, measurements, correlations, nuisance_correlations)

posterior = PosteriorDensity(model)

algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)
samples = bat_sample(posterior, algorithm).result;

p = plot(samples)
savefig(p, "plot.pdf")

# Warning:
# When using nuisance correlations, it may happen that the fit does not converge anymore.
# This happens mainly when the correlation values are close to -1 or +1.
# In such a case, you can try to reduce the allowed ranges in the priors for the Nuisance
# correlations to investigate at which values this happens.


# ============= Ranking of measurements and uncertainties ========================#
# With the `rank_measurements` and `rank_uncertainties` functions, the influence of
# the individual measurements or uncertainty types on the result of a fit can be estimated.
# For the ranking, each active measurement (respectively uncertainty type) is deactivated
# at a time and the fit is repeated. The results of all fits with one deactivated
# measurement (or uncertainty type) are then compared to the fit result with all measurements (uncertainty types) activated.
# A ranking is calculated based on a ranking criterion calculated from the posterior distributions of these fits.

# The default ranking criterion is the relative increase of the total width of the
# smallest interval containing 90% of the posterior probability when deactivating a measurement.
# For models with more than one parameter, the sum of relative increases of all
# one-dimensional smallest intervals is used. The default ranking criterion is therefore
#`SumOfSmallestIntervals(p=0.9, bins=200)`.
measurement_ranking = EFTfitter.rank_measurements(model)

# For ranking the uncertainty types, the relative decrease is used.
uncertainty_ranking = EFTfitter.rank_uncertainties(model, criterion = SumOfSmallestIntervals(p=0.9, bins=200))

# Please see the [ranking documentation] for further ranking criteria and keyword arguments.

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

