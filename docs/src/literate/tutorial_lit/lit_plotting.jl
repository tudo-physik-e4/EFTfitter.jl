#!jl # # EFTfitter.jl - Plotting Tutorial
#jl # EFTfitter.jl - Plotting Tutorial

#~EFTfitter includes several recipes for plotting its datatypes 
#~using [Plots.jl](http://docs.juliaplots.org/latest/)
using EFTfitter
using BAT
using IntervalSets
using Distributions
using Plots

#~we use the inputs from the basic tutorial:
include("tutorial_inputs.jl")
model = EFTfitterModel(parameters, measurements, correlations)

#~Note: All plots generated with the following plot recipes can be customized using the
#~[series attributes](http://docs.juliaplots.org/latest/generated/attributes_series/),
#~[axis attributes](http://docs.juliaplots.org/latest/generated/attributes_axis/),
#~[subplot attributes](http://docs.juliaplots.org/latest/generated/attributes_subplot/) and
#~[plot attributes](http://docs.juliaplots.org/latest/generated/attributes_plot/)


#!jl # ## Plotting Observables
#jl #~---------- Plotting Observables -------------------------

#~Plotting an `Observable` object:
plot(Observable(xsec1), (C1=0, C2=-1:0.01:1))

#~When plotting an `Observable` from the `EFTfitterModel`, it can be accessed in different ways:
plot(get_observables(model).xsec1, (C1=0, C2=-1:0.01:1))

plot(get_measurements(model).Meas1.observable, (C1=0, C2=-1:0.01:1))

#src savefig(plot(get_measurements(model).Meas1.observable, (C1=0, C2=-1:0.01:1)), "observable_plot.png")
#md # ![example plot](plots/observable_plot.png)

#~If the model has many parameters, it can be convenient to pass the paramter that should be 
#~plotted together with as a `NamedTuple` with default values for all parameters.
default_parameters = (C1=1, C2=0)
plot(get_observables(model).xsec1, (C2=-1:0.01:1,), default_parameters)
#src savefig(plot(get_observables(model).xsec1, (C2=-1:0.01:1,), default_parameters), "observable_plot_2.png")
#md # ![example plot](plots/observable_plot_2.png)

#~The second argument in this function overwrites the corresponding default parameters, 
#~so it is also possible to pass multiple parameters:
plot(get_observables(model).xsec1, (C2=-1:0.01:1, C1=2.3), default_parameters)
#src savefig(plot(get_observables(model).xsec1, (C2=-1:0.01:1, C1=2.3), default_parameters), "observable_plot_3.png")
#md # ![example plot](plots/observable_plot_3.png)

#~All observables of a model can easily be plotted in one plot:
p = plot()
for meas in get_measurements(model)
    p=plot!(meas.observable, (C1=0, C2=-1:0.01:1), ylabel="prediction")
end
p
#src savefig(p, "observable_plot_4.png")
#md # ![example plot](plots/observable_plot_4.png)

# When plotting observables, the default title contains the values of the fixed
# parameters. In case the title is too long for one line, linebreaks can be inserted
# using the keyword `titlewidth`. e.g.:
plot(get_observables(model).xsec1, (C1=-10:0.01:10, C2=0, C3=100, C4=200), titlewidth=13)
#src savefig(plot(get_observables(model).xsec1, (C1=-10:0.01:10, C2=0, C3=100, C4=200), titlewidth=13), "observable_plot_5.png")
#md # ![example plot](plots/observable_plot_5.png)

#!jl # ## Plotting Measurements
#jl #~---------- Plotting Measurements -------------------------

#~`Measurement` objects can be plotted on top of the observables as a horizontal line with an uncertainty band:
plot(get_measurements(model).Meas1.observable, (C1=0, C2=-0.2:0.01:0.2))
plot!(measurements.Meas1)
#src p = plot(get_measurements(model).Meas1.observable, (C1=0, C2=-0.2:0.01:0.2))
#src p = plot!(measurements.Meas1)
#src savefig(p, "measurement_plot_1.png")
#md # ![example plot](plots/measurement_plot_1.png)


#~However, when plotting the measurements of the `EFTfitterModel`, the following syntax 
#~is preferred as it supports showing the names of the measurments in the legend:
plot(get_measurements(model).Meas1.observable, (C1=0, C2=-0.2:0.01:0.2))
plot!(get_measurements(model), :Meas1)
#src p = plot(get_measurements(model).Meas1.observable, (C1=0, C2=-0.2:0.01:0.2))
#src p = plot!(get_measurements(model), :Meas1)
#src savefig(p, "measurement_plot_2.png")
#md # ![example plot](plots/measurement_plot_2.png)

#~The uncertainty typed to be plotted can be specified:
plot(get_measurements(model).Meas1.observable, (C1=0, C2=-0.2:0.01:0.2))
plot!(get_measurements(model), :Meas1, uncertainties=(:stat, :another_unc))
#src p = plot(get_measurements(model).Meas1.observable, (C1=0, C2=-0.2:0.01:0.2))
#src p = plot!(get_measurements(model), :Meas1, uncertainties=(:stat, :another_unc))
#src savefig(p, "measurement_plot_3.png")
#md # ![example plot](plots/measurement_plot_3.png)

#~When mutliple types of uncertainties are given, the sum of the squares is used as the total uncertainty.
#~By default, all uncertainties included in the `EFTfitterModel` are used.

#!jl # ## Plotting MeasurementDistributions
#jl #~---------- Plotting MeasurementDistributions -------------------------
#~`MeasurementDistribution`s can be plotted for fixed parameters:
plot(get_measurement_distributions(model).MeasDist.observable, (C1=1.2, C2=0))
plot!(get_measurement_distributions(model), :MeasDist)
#src p = plot(get_measurement_distributions(model).MeasDist.observable, (C1=1.2, C2=0))
#src p = plot!(get_measurement_distributions(model), :MeasDist)
#src savefig(p, "measdist_plot_1.png")
#md # ![example plot](plots/measdist_plot_1.png)

#~alternative plotting style for measurement distributions:
plot(get_measurement_distributions(model).MeasDist.observable, (C1=1.2, C2=0))
plot!(get_measurement_distributions(model), :MeasDist, st=:scatter)
#src p = plot(get_measurement_distributions(model).MeasDist.observable, (C1=1.2, C2=0))
#src p = plot!(get_measurement_distributions(model), :MeasDist, st=:scatter)
#src savefig(p, "measdist_plot_2.png")
#md # ![example plot](plots/measdist_plot_2.png)

#~Also for `MeasurementDistribution`s the uncertainty types to be plotted can be specified.
#~The names of the bins can be customized using the `bin_names` keyword.
plot(get_measurement_distributions(model).MeasDist.observable, (C1=1.2, C2=0))
plot!(get_measurement_distributions(model), :MeasDist, st=:scatter, uncertainties=(:stat,), bin_names=("First bin", "Second bin"))
#src p = plot(get_measurement_distributions(model).MeasDist.observable, (C1=1.2, C2=0))
#src p = plot!(get_measurement_distributions(model), :MeasDist, st=:scatter, uncertainties=(:stat,), bin_names=("First bin", "Second bin"))
#src savefig(p, "measdist_plot_3.png")
#md # ![example plot](plots/measdist_plot_3.png)
