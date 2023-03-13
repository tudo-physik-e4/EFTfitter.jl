using BAT
using EFTfitter
using Distributions
using DensityInterface
using Plots
gr(size=(1000,800), thickness_scaling=2)



# Counting experiment (Poisson distribution)
# Expected Background: b = 3.4
# Observed: N=3
# What is 95% upper limit on s?


likelihood = logfuncdensity(params -> begin
    s = params.s
    b = params.b
    n = 0#3

    return logpdf.(Poisson(s+b), n)#logpdf(Normal(1.9, 1.6), s)#logpdf.(Poisson(s+b), n)

end)

prior = BAT.NamedTupleDist(
    s = Uniform(0, 12),
    b = 0#3.4, #truncated(Normal(3.4, 2.), 0.01, Inf)
)

posterior = PosteriorMeasure(likelihood, prior);

samples, chains = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 4*10^6));

plot(samples, intervals=[0.9], bins=500)

interval_edges = get_smallest_interval_edges(samples, :s, 0.9, bins=500)

upper_limit = interval_edges.upper[1] 

using SpecialFunctions
function Normal_from_limit(best_fit_value, limit, confidence_level)
    μ = best_fit_value
    p = confidence_level
    q = limit

    σ = (q - μ)/(sqrt(2)*erfinv(p))

    return Normal(μ, σ)
end 

function Exponential_from_limit(limit, confidence_level)
    p = confidence_level
    q = limit

    λ = -1/q * log(1-p)

    return Exponential(1/λ)
end 

exponential_dist = Exponential_from_limit(upper_limit, 0.9)
normal_dist = Normal_from_limit(-10, upper_limit, 0.9)

x = 0:0.1:12
plot(samples, :s, intervals=[0.9], bins=500, label="Posterior")
plot!(x, pdf.(exponential_dist, x), lw=4, lc=1, label="Exponential from limit", legend=true)
plot!(x, 20*pdf.(normal_dist, x), lw=4, lc=2, label="Normal from limit", legend=true)


g(x) = (x+3.4)^3/6 * exp(-(x+3.4))

g(x) =   exp(-(x))
plot!(x, 1.8*g.(x), lw=4, lc="black", label="Poisson", legend=true)


g(x; μ=0, N=0) = (x+μ)^N * exp(-(x+μ))

x = 0:0.05:6
plot(x, g.(x, μ=2, N=4))

#------------------------------------------------------------
using DelimitedFiles

data_file = "C:\\Users\\Cornelius\\Projects\\plot2.csv"

data = readdlm(data_file, ';', Float64)

idxs = sortperm(data[:,1])
points_x = data[:, 1][idxs]
points_y = 1 .- data[:, 2][idxs]
points_y2 = data[:, 2][idxs]

x = 0.22:0.01:5.6
x1 = 0.22:0.01:2
x2 = 2.01:0.01:5.8

plot(points_x, points_y, st=:scatter, ms=2)

using Interpolations
itp = interpolate((points_x,), points_y, Gridded(Linear()))

plot!(x, itp.(x1), lw=2)


dx = only.(Interpolations.gradient.(Ref(itp), x))

plot!(x, dx, lw=2)

poi(N) = x-> x^N/factorial(N)*exp(-x)

plot(points_x, points_y2, st=:scatter, msw=0, ms=2)
#plot!(x1, 2*poi(0).(x1))
a=0.64
a2 = 2.2
plot!(x1, 0.35 .+ a.*exp.(-a*x1))
plot!(x2, -0.05 .+ a2.*exp.(-a*x2))

f1(x; a=0.64) = 1-(0.35 .+ a.*exp.(-a*x))
f2(x; a=2.2) = 1-(-0.05 .+ a.*exp.(-a*x))

plot(points_x, points_y, st=:scatter, ms=2, msw=0)
plot!(x1, f1.(x1), lw=2)
plot!(x2, f2.(x2), lw=2)



#---------------------------------------
n = 47
plot(points_x[1:n], points_y[1:n], st=:scatter, msw=0, ms=2)
plot!(points_x[n:end], points_y[n:end], st=:scatter, msw=0, ms=2, color=2)

@. ffq(x, p) = p[1]*x^2 + p[2]*x + p[3] + p[4]*x^3 + p[5]*x^4


p0 = ones(5)
fit = curve_fit(ffq, points_x[1:n], points_y[1:n], p0)
p1 = fit.param
plot!(x1, ffq(x1, p1), lw=3)


fit = curve_fit(ffq, points_x[n:end], points_y[n:end], p0)
p2 = fit.param
plot!(x2, ffq(x2, p2), lw=3)


ffq(p) = x -> ffq(x, p) 

using ForwardDiff
plot!(x1, ForwardDiff.derivative.(ffq(p1), x1))
plot!(x2, ForwardDiff.derivative.(ffq(p2), x2))

plot!(x2, pdf.(Normal(2, 1.6), x2))













using StatsFuns

plot(points_x, points_y, st=:scatter, msw=0, ms=2)


using LsqFit

@. model2(x, p) =  p[3]*StatsFuns.nchisqpdf(x, p[1], p[2])

p0 = [0.05, 1, 0.1]  # Initial guesses for the two parameters of the model
# Fitting the model to the data using nonlinear least-squares optimization
fit = curve_fit(model2, points_x[1:20], points_y[1:20], p0)
fit.param

plot!(x1, model2(x1, fit.param))