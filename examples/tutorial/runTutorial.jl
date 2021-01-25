# EFTfitter.jl - Tutorial
# This tutorial introduces the basic functionalities of EFTfitter.jl using a generic example.
# More functionalities of EFTfitter.jl, like handling nuisance correlations
# or ranking measurements and uncertainties, are shown in the
# [advanced tutorial](https://tudo-physik-e4.github.io/EFTfitter.jl/dev/advanced_tutorial/).

# Here, we create the `EFTfitterModel` from our inputs and run the actual analysis.

# First, we need to setup EFTfitter, BAT and some other Julia packages:
using EFTfitter
using BAT            # for sampling
using IntervalSets   # for specifying the prior
using Distributions  # for specifying the prior
using Plots # for plotting

# We include the definitions of observables, measurements,
# uncertainties and correlations from our `tutorial_inputs.jl` file:
include("tutorial_inputs.jl")

# We can then build the `EFTfitterModel` which combines all our inputs into
# one object that is then used to perform the analysis on.
model = EFTfitterModel(parameters, measurements, correlations)

# To sample the posterior distribution, we specify that our `EFTfitterModel`
# should be used and then setup BAT.jl to sample the EFTfitter likelihood.
posterior = PosteriorDensity(model)

algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)
samples = bat_sample(posterior, algorithm).result;
# For further information on settings & algorithms when sampling with BAT.jl
# see the BAT.jl [tutorial](https://bat.github.io/BAT.jl/dev/tutorial/#Parameter-Space-Exploration-via-MCMC)
# and [documentation](https://bat.github.io/BAT.jl/dev/stable_api/#BAT.bat_sample).

# We can then inspect the results of the sampling using BAT.jl's `SampledDensity`,
# giving a summary of the sampling and the results of the model parameters.
sd = SampledDensity(posterior, samples)
display(sd)

# Of course, plotting the resulting posterior distributions is also simple
# using Plots.jl and the BAT.jl plotting recipes:
p = plot(samples)
savefig(p, "plot.pdf")

# For information on how to customize plots of the samples, please see the BAT.jl
# [plotting documentation](https://bat.github.io/BAT.jl/dev/plotting/) and
# [examples](https://github.com/bat/BAT.jl/blob/master/examples/dev-internal/plotting_examples.jl).

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

