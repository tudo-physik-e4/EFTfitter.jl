#!jl # # EFTfitter.jl - Tutorial
#jl # EFTfitter.jl - Tutorial

#jl #~This tutorial introduces the basic functionalities of EFTfitter.jl using a generic example.  
#jl #~More functionalities of EFTfitter.jl, like handling nuisance correlations
#jl #~or ranking measurements and uncertainties, are shown in the  
#jl #~[advanced tutorial](https://tudo-physik-e4.github.io/EFTfitter.jl/dev/advanced_tutorial/). 
#md # Table of contents:
#md # ```@contents
#md # Pages = ["tutorial.md"]
#md # Depth = 3
#md # ```
#nb #~First, we need to setup EFTfitter, BAT and some other Julia packages:
#nb using EFTfitter
#nb using BAT           # for sampling
#nb using IntervalSets  # for specifying the prior
#nb using Distributions # for specifying the prior
#nb using Plots         # for plotting

#!jl include("lit_tutorial_inputs.jl")
#!jl # ## File "runTutorial.jl"
#~Here, we create the `EFTfitterModel` from our inputs and run the actual analysis.

#!nb #~First, we need to setup EFTfitter, BAT and some other Julia packages:
#!nb using EFTfitter
#!nb using BAT            # for sampling
#!nb using IntervalSets   # for specifying the prior
#!nb using Distributions  # for specifying the prior
#!nb using Plots # for plotting
#jl  
#jl #~We include the definitions of observables, measurements, 
#jl #~uncertainties and correlations from our `tutorial_inputs.jl` file:
#jl include("tutorial_inputs.jl")

#~We can then build the `EFTfitterModel` which combines all our inputs into 
#~one object that is then used to perform the analysis on.
model = EFTfitterModel(parameters, measurements, correlations)

#~To sample the posterior distribution, we specify that our `EFTfitterModel` 
#~should be used and then setup BAT.jl to sample the EFTfitter likelihood. 
posterior = PosteriorMeasure(model)

algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)
samples = bat_sample(posterior, algorithm).result;
#~For further information on settings & algorithms when sampling with BAT.jl 
#~see the BAT.jl [tutorial](https://bat.github.io/BAT.jl/dev/tutorial/#Parameter-Space-Exploration-via-MCMC) 
#~and [documentation](https://bat.github.io/BAT.jl/dev/stable_api/#BAT.bat_sample).

#~We can then inspect the results of the sampling using BAT.jl's `bat_report`, 
#~giving a summary of the sampling and the results of the model parameters.
bat_report(samples)
#md ```
#md ```
#md BAT.jl - SampledDensity
#md ──────────────────────────────
#md 
#md Sampling:
#md ─────────────────────────
#md total number of samples:      103343
#md effective number of samples: (C1 = 19897.99484950266, C2 = 18106.884734831812)
#md 
#md 
#md Parameter estimates:
#md ─────────────────────────
#md number of free parameters: 2
#md 
#md Table with 5 columns and 2 rows:
#md      parameter  mean       std        global_mode  marginal_mode
#md    ┌────────────────────────────────────────────────────────────
#md  1 │ C1         0.864514   0.349673   0.821192     0.81
#md  2 │ C2         0.0129843  0.0210368  0.0155112    0.0165
#md 
#md 
#md Covariance matrix:
#md ─────────────────────────
#md 2×2 Named Array{Float64,2}
#md cov ╲  │          C1           C2
#md ───────┼─────────────────────────
#md C1     │    0.122271   -0.0070394
#md C2     │  -0.0070394  0.000442548

#~Information about the smallest 1d intervals containing p% proability can be
#~obtained using the `get_smallest_interval_edges` function:
intervals_1d_C1 = get_smallest_interval_edges(samples, :C1, 0.9, bins=200, atol=0.1)
println("lower interval edges: $(intervals_1d_C1.lower)")
println("upper interval edges: $(intervals_1d_C1.upper)")

#~The keyword `atol` controls the absolute tolerance for which intervals are joined 
#~together when they are seperated less than this value. This is particularly useful
#~when a large number of bins is used.

#~Of course, plotting the resulting posterior distributions is also simple 
#~using Plots.jl and the BAT.jl plotting recipes:
p = plot(samples)
#!nb savefig(p, "plot.pdf")
#md # ![example plot](plots/plot.png)

#~For information on how to customize plots of the samples, please see the BAT.jl 
#~[plotting documentation](https://bat.github.io/BAT.jl/dev/plotting/) and 
#~[examples](https://github.com/bat/BAT.jl/blob/master/examples/dev-internal/plotting_examples.jl).

p = plot(samples, 0.9)
#!nb savefig(p, "plot_1d.pdf")
#md # ![example plot](plots/interval_plot_1.png)

#~For customizing the plots of the 1D intervals, also see the EFTfitter 
#~[plotting documentation](https://tudo-physik-e4.github.io/EFTfitter.jl/dev/plotting/) 
#~and [tutorial](https://github.com/tudo-physik-e4/EFTfitter.jl/blob/main/examples/tutorial/plotting.jl).
