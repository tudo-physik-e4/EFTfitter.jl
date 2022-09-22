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
posterior = PosteriorMeasure(model)

algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)
samples = bat_sample(posterior, algorithm).result;
# For further information on settings & algorithms when sampling with BAT.jl
# see the BAT.jl [tutorial](https://bat.github.io/BAT.jl/dev/tutorial/#Parameter-Space-Exploration-via-MCMC)
# and [documentation](https://bat.github.io/BAT.jl/dev/stable_api/#BAT.bat_sample).

# We can then inspect the results of the sampling using BAT.jl's `SampledDensity`,
# giving a summary of the sampling and the results of the model parameters.
sd = SampledDensity(posterior, samples)
display(sd)

# Information about the smallest 1d intervals containing p% proability can be
# obtained using the `get_smallest_interval_edges` function:
intervals_1d_C1 = get_smallest_interval_edges(samples, :C1, 0.9, bins=200, atol=0.1)
println("lower interval edges: $(intervals_1d_C1.lower)")
println("upper interval edges: $(intervals_1d_C1.upper)")

# The keyword `atol` controls the absolute tolerance for which intervals are joined
# together when they are seperated less than this value. This is particularly useful
# when a large number of bins is used.

# Of course, plotting the resulting posterior distributions is also simple
# using Plots.jl and the BAT.jl plotting recipes:
p = plot(samples)
savefig(p, "plot.pdf")

# For information on how to customize plots of the samples, please see the BAT.jl
# [plotting documentation](https://bat.github.io/BAT.jl/dev/plotting/) and
# [examples](https://github.com/bat/BAT.jl/blob/master/examples/dev-internal/plotting_examples.jl).

p = plot(samples, 0.9)
savefig(p, "plot_1d.pdf")

# For customizing the plots of the 1D intervals, also see the EFTfitter
# [plotting documentation](https://tudo-physik-e4.github.io/EFTfitter.jl/dev/plotting/)
# and [tutorial](https://github.com/tudo-physik-e4/EFTfitter.jl/blob/main/examples/tutorial/plotting.jl).

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

