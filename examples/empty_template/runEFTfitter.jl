# EFTfitter.jl - Empty Template

using EFTfitter
using BAT            # for sampling
using IntervalSets   # for specifying the prior
using Distributions  # for specifying the prior
using Plots # for plotting

# include definitions of observables, measurements,
# uncertainties and correlations from the inputs file:
include("inputs.jl")

# create an `EFTfitterModel` object:
model = EFTfitterModel(parameters, measurements, correlations)
# create posterior distribution:
posterior = PosteriorMeasure(model);

# sample the posterior distribution with BAT.jl:
algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)
samples = bat_sample(posterior, algorithm).result;

# let's get a quick overview of results:
bat_report(samples)

#  plot the posterior distribution:
p = plot(samples)
savefig(p, "plot.pdf")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
