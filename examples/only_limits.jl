using EFTfitter
using BAT            
using IntervalSets   
using Distributions  
using Plots 

parameters_BR = BAT.NamedTupleDist(
    BR_tcH = 0.0..3, 
    BR_tuH = 0.0..3,
)

parameters_EFT = BAT.NamedTupleDist(
    Ccϕ = 0.0..2.0, 
    Cuϕ = 0.0..2.0,
)

fakefunc(params) = 1000.

BR_cH_EFT(params) = params.Ccϕ^2 / 1946.6
BR_uH_EFT(params) = params.Cuϕ^2 / 1946.6

BR_cH(params) = params.BR_tcH
BR_uH(params) = params.BR_tuH

measurements = (
    meas1 = Measurement(fakefunc, 123., uncertainties = (unc1=0.1, ), active=true),
    #Hu = Measurement(Observable(BR_uH_EFT, min=0.), 0.37e-3, uncertainties = (unc1=0.18e-3, ), active=false),
    #Hc = Measurement(Observable(BR_cH_EFT, min=0.), 0.51e-3, uncertainties = (unc1=0.24e-3, ), active=false),
)


limits_EFT = (
    Hu = ExponentialUpperLimit(BR_uH_EFT, 6.9e-4, 0.95, active=true),
    Hc = ExponentialUpperLimit(BR_cH_EFT, 9.4e-4, 0.95, active=true),
)

limits_EFT_expected = (
    Hu = ExponentialUpperLimit(BR_uH_EFT, 3.5e-4, 0.95, active=true),
    Hc = ExponentialUpperLimit(BR_cH_EFT, 4.8e-4, 0.95, active=true),
)

limits_BR = (
    Hu = ExponentialUpperLimit(BR_uH, 6.9e-1, 0.95, active=true),
    Hc = ExponentialUpperLimit(BR_cH, 9.4e-1, 0.95, active=true),
)

limits_BR_expected = (
    Hu = ExponentialUpperLimit(BR_uH, 3.5e-1, 0.95, active=true),
    Hc = ExponentialUpperLimit(BR_cH, 4.8e-1, 0.95, active=true),
)


correlations = (
    unc1 = NoCorrelation(active=true),
)

model_EFT = EFTfitterModel(parameters_EFT, measurements, correlations, limits=limits_EFT)
model_EFT_expected = EFTfitterModel(parameters_EFT, measurements, correlations, limits=limits_EFT_expected)

model_BR = EFTfitterModel(parameters_BR, measurements, correlations, limits=limits_BR)
model_BR_expected = EFTfitterModel(parameters_BR, measurements, correlations, limits=limits_BR_expected)

posterior_EFT = PosteriorMeasure(model_EFT)
posterior_EFT_expected = PosteriorMeasure(model_EFT_expected)
posterior_BR = PosteriorMeasure(model_BR)
posterior_BR_expected = PosteriorMeasure(model_BR_expected)

algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 1*10^5, nchains = 4)

samples_EFT = bat_sample(posterior_EFT, algorithm).result;
samples_EFT_expected = bat_sample(posterior_EFT_expected, algorithm).result;
samples_BR = bat_sample(posterior_BR, algorithm).result;
samples_BR_expected = bat_sample(posterior_BR_expected, algorithm).result;


plot(samples_BR, left_margin=14Plots.px)




# Load ATLAS DATA
using CSV, DataFrames
atlas_EFT = CSV.File("..\\..\\atlas_EFT.csv", header=false) |> DataFrame
atlas_EFT_expected = CSV.File("..\\..\\atlas_EFT_expected.csv", header=false) |> DataFrame
atlas_BR = CSV.File("..\\..\\atlas_BR.csv", header=false) |> DataFrame
atlas_BR_expected = CSV.File("..\\..\\atlas_BR_expected.csv", header=false) |> DataFrame

#----- Plot BR ----------------------------------------------------------------------------------
f=100
# ATLAS Data
p = plot(atlas_BR.Column1*1e3, atlas_BR.Column2*1e3, lw=15, label="ATLAS observed limit", lc=:black, alpha=0.4)
p = plot!(atlas_BR_expected.Column1*1e3, atlas_BR_expected.Column2*1e3, lw=15, label="ATLAS expected limit", c=:black, st=:scatter, alpha=0.4)

# Samples
#p = plot!(samples_BR, (1,2), st=:smallest_intervals, intervals=[0.95,], bins=100, color=:red, alpha=0.1, marginalmode=false)

# EFTfitter bounds
p = plot!(samples_BR, (1,2), st=:smallest_intervals_contour, intervals=[0.95,], bins=20, smoothing=0.1, lw=5, lc=:red, 
marginalmode=false,  dpi=600,  size=(f*8.3, f*5.7), tickfontsize=14,  label="EFTfitter")

p = plot!(samples_BR_expected, (1,2), st=:smallest_intervals_contour, intervals=[0.95,], bins=20, smoothing=0.1, lw=5, alpha=0.5, lc=:red, 
marginalmode=false,  dpi=600, size=(f*8.3, f*5.7), tickfontsize=14,  label="EFTfitter")

# Label
#p = plot!([0], [-10], lw=3, label="EFTfitter", lc=:red, xlims=(0, 1.2e-3), ylims=(0, 1.2e-3),
#    xticks=0:0.2e-3:1.2e-3, ytick=0:0.2e-3:1.2e-3, xlabel="BR(t->cH)", ylabel="BR(t->uH)")

savefig(p, "..\\..\\exponential_BR.png")

p



#----- Plot EFT ----------------------------------------------------------------------------------
f=100
# ATLAS Data
p = plot(atlas_EFT.Column1, atlas_EFT.Column2, lw=15, label="ATLAS observed limit", lc=:black, alpha=0.4)
p = plot!(atlas_EFT_expected.Column1, atlas_EFT_expected.Column2, lw=15, label="ATLAS expected limit", c=:black, st=:scatter, alpha=0.4)

# Samples
#p = plot!(samples_EFT, (1,2), st=:smallest_intervals, intervals=[0.95,], bins=100, color=:red, alpha=0.1, marginalmode=false)

# EFTfitter bounds
p = plot!(samples_EFT, (1,2), st=:smallest_intervals_contour, intervals=[0.95,], bins=20, smoothing=0.1, lw=5, lc=:red, 
marginalmode=false,  dpi=600, xticks=0:0.2:1.8, ytick=0:0.2:2.0, size=(f*8.3, f*5.7), tickfontsize=14,  label="EFTfitter")

p = plot!(samples_EFT_expected, (1,2), st=:smallest_intervals_contour, intervals=[0.95,], bins=20, smoothing=0.1, lw=5, alpha=0.5, lc=:red, 
marginalmode=false,  dpi=600, xticks=0:0.2:1.8, ytick=0:0.2:2.0, size=(f*8.3, f*5.7), tickfontsize=14,  label="EFTfitter")

# Label
p = plot!([0], [-10], lw=3, label="EFTfitter", lc=:red, xlims=(0, 1.8), ylims=(0, 2))

savefig(p, "..\\..\\exponential_EFT.png")

p


#---------------- Test -------------------------
λ1 = - log(1. - 0.95)/9.4e-4
d1 = Exponential(1/λ1)

λ2 = - log(1. - 0.95)/6.9e-4
d2 = Exponential(1/λ2)


















marg = BAT.MarginalDist(samples, :Cuϕ, bins=500)
marghist = convert(Histogram, marg.dist)


using FHist
fh = Hist1D(marghist)

integral(fh)

sum = 0. 
inte = integral(fh)
for b in 1:nbins(fh)
    sum += bincounts(fh)[b]
    if sum/inte > 0.95
        println("lower: ", bincenters(fh)[b])
        break
    end
end






using StatsBase

function get_smallest_interval_edges(
    samples::DensitySampleVector,
    key::Union{Symbol,Real},
    p::Real;
    bins=200,
    atol=0.0)

    marg = BAT.MarginalDist(samples, key, bins=bins)

    marghist = convert(Histogram, marg.dist)
    hist_p = BAT.get_smallest_intervals(marghist, [p])[1][1]

    lower, upper = BAT.get_interval_edges(hist_p, atol=atol)
    return (lower=lower, upper=upper)
end

_lower, _upper = get_smallest_interval_edges(samples_BR, :BR_tuH, 0.95, bins=400, atol=0.01)
_lower, _upper = get_smallest_interval_edges(samples_BR, :BR_tcH, 0.95, bins=400, atol=0.01)


plot(samples, :Ccϕ, bins=400)

d1 = posterior.likelihood.density._d.limit_distributions[1]

x1 = 0.1e-6:0.1e-5:7e-4
plot(x1, pdf.(d1, x1), label="limit distribution")


plot(rand(d1, 10000), st=:hist)

lim = 6.9e-4
p = 0.95

λ = -log(1-p) / lim

d = Exponential(1/λ)

cdf(d, 6.9e-4)
quantile(d, 0.95)

Z = 1.96 # this is for 95% CL
μ = 0
X = 6.9e-4
σ = (X-μ)/Z 

d = Normal(μ, σ)
quantile(truncated(d, 0, Inf), 0.95)