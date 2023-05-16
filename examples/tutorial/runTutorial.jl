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
get_total_covariance(model)
run_speed_test(model)

# To sample the posterior distribution, we specify that our `EFTfitterModel`
# should be used and then setup BAT.jl to sample the EFTfitter likelihood.
posterior = PosteriorMeasure(model)
posterior.likelihood.density._d.observable_mins

using DensityInterface
v=rand(parameters)
@btime logdensityof(posterior)(v)
algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)
samples = bat_sample(posterior, algorithm).result;

samples.v
@btime samples.v.C1
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

struct ResultT
    a::Float64
    b::Float64
end

ResultT(a::Float64) = ResultT(a, 0.)

ResultT(a::Tuple) = ResultT(a[1], a[2])

function myf()
    return 1.
end 
function myf(a)
    return 1.
end 

function myf2()
    return 1., 2.
end

function myf2(a)
    return 1., 2.
end

function myf3()
    return ResultT(1., 2.)
end

function evalfs(fs)
    return [f(3) for f in fs]
end
        
preds = zeros(5)
uncs = zeros(5)

function evalfs2(fs, preds=preds, uncs=uncs)
    for i in eachindex(fs)
        p, u = fs[i]()
        preds[i] = p
        uncs[i] = u
    end
    return preds, uncs
end


function evalfs2e(fs, preds::Vector{Float64}=preds, uncs::Vector{Float64}=uncs)
    for i in eachindex(fs)
        #println(fs[i]())
        r = ResultT((1.,2.))#fs[i]())
        #println(r)
        preds[i] = r.a
        uncs[i] = r.b
    end

    return preds, uncs
end

function evalfs3(fs, preds=preds, uncs=uncs)
    for i in eachindex(fs)
        r = fs[i](5)
        preds[i] = r[1]
        uncs[i] = r[2]
    end
    return preds, uncs
end

function ff(f::Function)
    return x->(f(x), 0.)
end


g = x -> (myf(x), 0.)

g(5)


@btime a = evalfs([myf, myf, myf, myf2, myf2])

@btime a = evalfs3([ff(myf), ff(myf), myf2, myf2, myf2])

@btime c = evalfs2([myf2, myf2, myf2, myf2, myf2])

@btime c = evalfs3([myf3, myf3, myf3, myf3, myf3])


@profview [evalfs2b([myf2, myf2, myf2, myf2, myf2]) for i in 1:100]

@code_warntype evalfs2e([myf2, myf2, myf2, myf2, myf])

@code_warntype evalfs3([ff(myf), ff(myf), myf2, myf2, myf2])


using DensityInterface
logdensityof(posterior)(rand(parameters))


using BenchmarkTools
using SparseArrays, LinearAlgebra
using BAT.TypedTables
using EFTfitter.Setfield


tbl, bms = run_speed_test(model, verbose=true)

run_speed_test(model, verbose=true)




bms[3]

tbl.Type[1]

vs3 = 0# rand(model.parameters, 100)
@btime(run_btime(posterior, vs3), setup=(vs3 = rand(model.parameters, 100)))


a = @benchmark rand(model.parameters)
median(a)
median(a.times)
a.memory
a.allocs
dump(a)

@btime(sort!(x), setup=(x = rand(5)), evals=1)

typeof(posterior.likelihood.density._d.invcov)

false ? println("Hi") : nothing



rarr2 = [[rand(), rand()] for i in 1:n]


function f1(n)
    r = []
    for i in 1:n
        if rand() > 0.5
            push!(r, (rand(), rand()))
        else
            push!(r, rand())
        end
    end
    return r
end

function f1nt(n)
    r = []
    for i in 1:n
        if rand() > 0.5
            push!(r, (prediction=rand(), unc=rand()))
        else
            push!(r, (prediction = rand(),))
        end
    end
    return r
end

function f2(n)
    r = [rand() for i in 1:n]
    r2 = [rand() for i in 1:n]
    return r, r2
end

n =  50000
r = f1(n)
rnt = f1nt(n)
r1, r2 = f2(n)

rarr = [[ri...] for ri in r]

function fs(a)
    return a[1]
end
using BenchmarkTools
@btime rs = [fs(ri) for ri in r]

@btime r = [ri.prediction for ri in rnt]

r

rnt.prediction

for i in 1:10000
    (p->p.prediction).(rnt)
end
@btime map(p->p.prediction, rnt)
@btime getfield.(rnt, :prediction)

@btime map(p->p[1], rarr)

@btime rarr[1:end]

z = @view map(p->p.prediction, rnt)

z = view(r[:], 1)


using ValueShapes
using BAT.ArraysOfArrays
rnt

@btime shape = valshape(rnt[1])
@btime data = nestedview(Array{Float64}(undef, totalndof(shape), length(rnt)))
@btime Y = shape.(data)
Y[:] = rnt

sna = ShapedAsNTArray(rarr2, shape)
sna.unc
@btime zeros(50000)



z = () -> 5

z()