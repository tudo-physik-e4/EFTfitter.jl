#!jl # # EFTfitter.jl - Empty Template
#jl # EFTfitter.jl - Empty Template

#nb using EFTfitter
#nb using BAT            # for sampling
#nb using IntervalSets   # for specifying the prior
#nb using Distributions  # for specifying the prior
#nb using Plots          # for plotting

#!jl include("lit_inputs.jl")

#!nb using EFTfitter
#!nb using BAT            # for sampling
#!nb using IntervalSets   # for specifying the prior
#!nb using Distributions  # for specifying the prior
#!nb using Plots # for plotting
#jl  
#jl #~include definitions of observables, measurements, 
#jl #~uncertainties and correlations from the inputs file:
#jl include("inputs.jl")

#~create an `EFTfitterModel` object:
model = EFTfitterModel(parameters, measurements, correlations)
#~create posterior distribution:
posterior = PosteriorDensity(model);

#~sample the posterior distribution with BAT.jl:
algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)
samples = bat_sample(posterior, algorithm).result;

#~create and display a `SampledDensity` object for a quick overview of results:
sd = SampledDensity(posterior, samples)
display(sd)

#~ plot the posterior distribution:
p = plot(samples)
#!nb savefig(p, "plot.pdf")
