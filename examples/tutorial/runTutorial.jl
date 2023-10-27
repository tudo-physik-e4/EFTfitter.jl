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
posterior = PosteriorMeasure(model)
#run_speed_test(model)

v = rand(parameters)
v = (C1 = 1.3646163105428428, C2 = 0.0669263861339656)

#@code_warntype EFTfitter.evaluate_funcs!(m.mus, m.observable_functions, v, m)

using BenchmarkTools
using DensityInterface
# logdensityof(posterior)(v)
@btime logdensityof(posterior)(v) # -71.62957930828858
# 350.000 ns (18 allocations: 1.19 KiB)



# To sample the posterior distribution, we specify that our `EFTfitterModel`
# should be used and then setup BAT.jl to sample the EFTfitter likelihood.

parameters1 = BAT.NamedTupleDist(
    p1 = -20..20, 
    p2 = -10..10,
)

function testfunc1(params)
    c = [20.12, 5.56, 325.556]
    return c[1] * params.p1^2 + c[2] * params.p1 * params.p2 + c[3] * params.p2^2
end

measurements1 = (
    meas1 = Measurement(testfunc1, 111.1, 
    uncertainties = (unc1=10.1, unc2=12.2, unc3=13.3), active=true),
    
    meas2 = Measurement(Observable(testfunc1, min=0, max=1000), 222.2, 
    uncertainties = (unc1=20.1, unc2=20.2, unc3=23.3), active=false),
    
    meas3 = Measurement(Observable(testfunc1, min=0, max=1000), 333.3, 
    uncertainties = (unc1=30.1, unc2=30.2, unc3=30.3), active=true),
    
    meas4 = MeasurementDistribution(Function[testfunc1, testfunc1, testfunc1],
    [10, 20, 30], uncertainties = (unc1=[0.11, 0.12, 0.13], unc2=[0.21, 0.22, 0.23], unc3=[0.31, 0.32, 0.33]),
    active = [true, false, true], bin_names=[Symbol("0_5"), Symbol("5_10"), Symbol("10_20")])
)

corr_matrix = to_correlation_matrix(measurements1,
    (:meas1, :meas2, 0.12), # will be ignored later in EFTfitterModel
    (:meas1, :meas3, 0.13), 
    (:meas1, :meas4_0_5, 0.141), 
    (:meas1, :meas4_5_10, 0.142), # will be ignored later in EFTfitterModel
    (:meas1, :meas4_10_20, 0.143), # will be ignored later in EFTfitterModel
    (:meas4_0_5, :meas4_5_10, 0.412), # will be ignored later in EFTfitterModel
    (:meas4_0_5, :meas4_5_10, 0.413), 
    (:meas4_0_5, :meas4_0_5, 0.9) # will be ignored later in EFTfitterModel
)

correlations1 = (
    unc1 = NoCorrelation(active=true),

    # wrong matrix size for number of measurements, will be ignored if active=false:
    unc2 = Correlation([1.0 0.5 0.7; 
                        0.5 1.0 0.6;
                        0.7 0.6 1.0], active=false),
    
    unc3 = Correlation(corr_matrix)
)

nuisance_correlations = (
    ρ1 = NuisanceCorrelation(:unc1, :meas1, :meas3, 0..0.5),
    ρ2 = NuisanceCorrelation(:unc1, :meas1, :meas2, truncated(Normal(0, 1), 0, 0.9)),
)

model1 = EFTfitterModel(parameters1, measurements1, correlations1, limits=nothing, nuisances=nuisance_correlations, CovarianceType=Matrix)





posterior = PosteriorMeasure(model1)
m = posterior.likelihood.density._d

v2 = (p1 = 1.3646163105428428, p2 = 0.0669263861339656, ρ1=0.2, ρ2=0.3 )

EFTfitter.get_current_invcov(m, v2)


#@btime EFTfitter.evaluate_funcs!(m.mus, m.observable_functions, v, m)
#187.771 ns (14 allocations: 800 bytes)




using DensityInterface
using BenchmarkTools
v = rand(parameters)

v = (C1 = 1.3646163105428428, C2 = 0.0669263861339656)

#@code_warntype EFTfitter.evaluate_funcs!(m.mus, m.observable_functions, v, m)

# logdensityof(posterior)(v)
@btime logdensityof(posterior)(v) # -71.62957930828858
# 350.000 ns (18 allocations: 1.19 KiB)
# 633.523 ns (18 allocations: 1.19 KiB)
# 983.333 ns (18 allocations: 1.19 KiB)

@btime [m.observable_functions[i](v) for i in eachindex(m.observable_functions)]

@code_warntype logdensityof(posterior)(v)
@code_llvm logdensityof(posterior)(v)

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
#--------------------------------s--------------



M = rand(4, 100)

m = view(M, 1, :)
x = rand(100)
r = m - x