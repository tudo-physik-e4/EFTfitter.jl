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

using SparseArrays
using LinearAlgebra
f(x) = Symmetric(sparse(x))

model = EFTfitterModel(parameters, measurements, correlations, CovarianceType=f)
#run_speed_test(model)

# To sample the posterior distribution, we specify that our `EFTfitterModel`
# should be used and then setup BAT.jl to sample the EFTfitter likelihood.
posterior = PosteriorMeasure(model)

using DensityInterface
using BenchmarkTools
v = rand(parameters)

v = (C1 = 1.3646163105428428, C2 = 0.0669263861339656)

@btime logdensityof(posterior)(v) # -71.62957930828858
# 983.333 ns (18 allocations: 1.19 KiB)

@code_warntype logdensityof(posterior)(v)

v = (C1 = 1.3777803296719995, C2 = -0.13471933089976204)
logdensityof(posterior)(v) # -139.56233066526895

v = (C1 = 0.4310417711590908, C2 = 0.5097850986277717)
logdensityof(posterior)(v) # -2162.171314291685


import Random
Random.seed!(1234)
for v in [rand(parameters) for i in 1:10000]
    println(logdensityof(posterior)(v))
end



@btime logdensityof(posterior)(v) #btime: 593.220 ns (26 allocations: 1.34 KiB)

algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4, strict=false)
include("ram_sampler.jl")
algorithm = RAMSampler(nchains=1, nsteps=6*10^5, nburnin=4*10^5)

@time samples = bat_sample(posterior, algorithm).result;


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

r = rand(500)
using LinearAlgebra
using BenchmarkTools
@btime Diagonal(r.^2)

N = 50
r = rand(N)*10
M = rand(N, N)

function addunc!(M, r)
    for i in 1:size(M)[1]
        M[i, i] += r[i]^2
    end
end

function addunc2(M, r)
   return M+Diagonal(r.^2)
end

@btime addunc!(M, r)

@btime addunc2(M, r)


M = [1 2; 3 4]
M2 = [1.5 2; 3 5]

invM = inv(M)
invM2 = inv(M2)

#------ Dispatch on Val -----------------------
f(a::Val{false}) = println("f=False")
f(a::Val{true}) = println("f=True")
f(a) = println("??")

@btime f(Val(false))
@btime f(3)



#----------------------------------------------
#- Testing Dispatch on Type of Field in struct

struct MyStruct{AT, BT}
    A::AT
    B::BT
    C::Float64
end

f(ms::MyStruct) = 1

f(ms::MyStruct{String, Float64}) = "Hi"

f(ms::MyStruct{Float64, String}) = "Ho"


m = MyStruct(1., 2., 3.)
m = MyStruct("H", 2., 3.)
m = MyStruct(3, "H0", 3.)

f(m)
#----------------------------------------------


