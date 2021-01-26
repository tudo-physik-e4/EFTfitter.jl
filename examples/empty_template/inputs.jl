# EFTfitter.jl - Empty Template
# ============= Parameters =============================================#
parameters = BAT.NamedTupleDist(
    p1 = -2..2,
)


# ============= Observables ============================================#
function observable1(params)
    return params.p1
end


# ============= Measurements ===========================================#
measurements = (
    Meas1 = Measurement(observable1, 0.0, uncertainties = (unc1 = 0.1,), active=true),

    #MeasDist = MeasurementDistribution(obs_array, values_array, uncertainties = (unc1 = unc1_array,), active=false),
)


# ============= Correlations ===========================================#
correlations = (
    unc1 = NoCorrelation(active=true),
)

#corr_matrix = to_correlation_matrix(measurements,
#  (:Meas1, :Meas2, 0.1),
#)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

