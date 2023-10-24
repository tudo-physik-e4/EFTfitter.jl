# EFTfitter.jl - Empty Template
using EFTfitter
using BAT            # for sampling
using IntervalSets   # for specifying the prior
using Distributions  # for specifying the prior
using Plots # for plotting

# include definitions of observables, measurements,
# uncertainties and correlations from the inputs file:
include("inputs.jl")
include("../tutorial/ram_sampler.jl")

# create an `EFTfitterModel` object:
model = EFTfitterModel(parameters, measurements, correlations)
# create posterior distribution:
posterior = PosteriorMeasure(model);
get_total_covariance(model)


# sample the posterior distribution with BAT.jl:
algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4, strict=false)
algorithm = RAMSampler(nchains=4, nsteps=8*10^5, nburnin=4*10^5)

samples = bat_sample(posterior, algorithm).result;

findmode = bat_findmode(posterior)
findmode.result

findmode.info
# create and display a `SampledDensity` object for a quick overview of results:
sd = SampledDensity(posterior, samples)
print(sd)
bat_report(sd)

cov(samples)

#  plot the posterior distribution:
p = plot(samples)
savefig(p, "plot.pdf")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

sqrt(180)

d = 6
0.5*(2*pi)^(d/2)
