# EFTfitter.jl - Curve Fit Example

using Random

# generate mock data
Random.seed!(45)
N = 20

f(x) = 1.1 + 0.2*x + 4.3*x^2 + 0.4*x^3 + 0.2*x^4 + 0.6*x^5
x=-10:0.1:10

x_data = rand(x, N)
ys = f.(x_data) 
y_data = [y + rand(Normal(0, 15)) for y in ys]



# ============= Parameters =============================================#
parameters = BAT.NamedTupleDist(
    p = [-20..40, -20..10, -10..10, -10..10, -5..5, -2..2],
    #p = fill(-1..1, 6),
)


# ============= Observables ============================================#
function observable(params, x)
    return params.p[1] + params.p[2]*x + params.p[3]*x^2 + params.p[4]*x^3 + params.p[5]*x^4 + params.p[6]*x^5
end

g(x) = params -> observable(params, x)

obs_array = Function[g(x) for x in x_data]

#unc_array = fill(1e-5, length(y_data))
unc_array = sqrt.(abs.(y_data))

# ============= Measurements ===========================================#
using LinearAlgebra
cov_matrix = Matrix(I, N, N)

cor, unc = cov_to_cor(cov_matrix)

unc = unc .* sqrt(180)


measurements = (
    MeasDist = MeasurementDistribution(obs_array, y_data, uncertainties = (unc1 = unc,)),
)


# ============= Correlations ===========================================#
correlations = (
    unc1 = NoCorrelation(active=true),
)

#corr_matrix = to_correlation_matrix(measurements,
#  (:Meas1, :Meas2, 0.1),
#)


# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

