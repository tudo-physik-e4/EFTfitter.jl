#jl # EFTfitter.jl - Empty Template 

#!jl # ### Parameters 
#jl #~============= Parameters =============================================#   
parameters = BAT.NamedTupleDist(
    p1 = -2..2, 
)


#!jl # ### Observables
#jl #~============= Observables ============================================#
function observable1(params)
    return params.p1
end


#!jl # ### Measurements
#jl #~============= Measurements ===========================================# 
measurements = (
    Meas1 = Measurement(observable1, 0.0, uncertainties = (unc1 = 0.1,), active=true),
    
    #MeasDist = MeasurementDistribution(obs_array, values_array, uncertainties = (unc1 = unc1_array,), active=false),  
)


#!jl # ### Correlations
#jl #~============= Correlations ===========================================#   
correlations = (
    unc1 = NoCorrelation(active=true),
)

#corr_matrix = to_correlation_matrix(measurements,
#    (:Meas1, :Meas2, 0.1), 
#)
