using Random
using Plots

# generate mock data
Random.seed!(45)
N = 20

f(x) = 1.1 + 0.2*x + 4.3*x^2 + 0.4*x^3 + 0.2*x^4 + 0.6*x^5
x=-10:0.1:10

x_data = rand(x, N)
ys = f.(x_data) 
y_data = [y + rand(Normal(0, 15)) for y in ys]

for xi in x_data
    print(xi, ",")
end

for yi in y_data
    print(yi, ",")
end



plot(x, f.(x), label="Truth")
plot!(x_data, y_data, st=:scatter, label="Data")


# using the Julia package "LsqFit.jl"
using LsqFit
@. fit_func(x, p) = p[1] + p[2]*x + p[3]*x^2 + p[4]*x^3 + p[5]*x^4 + p[6]*x^5

p0 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]  # initial guess

# Fitting the model to the data using nonlinear least-squares optimization
fit = curve_fit(fit_func, x_data, y_data, p0)




best_fit = fit.param  # Retrieving the best-fit parameter values for the model
sigma = stderror(fit)  # Calculating the standard error of the fit
margin_of_error = margin_error(fit, 0.05)  # Calculating the margin of error for the fit at a 95% confidence level
confidence_inter = confidence_interval(fit, 0.05)  # Calculating the confidence interval for the fit at a 95% confidence level

plot(x, f.(x), label="Truth")
plot!(x_data, y_data, st=:scatter, label="Data")
plot!(x, fit_func(x, best_fit), label="Fit")

sqrt(20)^