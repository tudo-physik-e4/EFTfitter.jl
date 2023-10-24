import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt


x = np.array([-6.5,-9.3,7.9,-8.9,9.7,5.5,5.5,2.5,-8.7,1.0,1.4,1.9,7.6,3.8,-8.5,-7.0,-2.8,0.9,1.3,4.2,])
y = np.array([437.903955201858,1537.749199238419,1231.9258895641094,1323.9224677984573,2527.6637713963096,399.57818043731953,362.4963567699898,56.85153464909146,1202.774601731368,25.205508693430282,17.192572423872164,1.914974838495624,1086.4978390611982,134.2200703399384,1107.3595441415214,561.0614957587213,49.062840977188266,20.906407232145973,9.70889789068209,144.96530373831723])

# Test function with coefficients as parameters
def test(x, a, b, c, d, e):
	return a + b*x + c*x**2 + d*x**3 + e*x**4

# curve_fit() function takes the test-function
# x-data and y-data as argument and returns
# the coefficients a and b in param and
# the estimated covariance of param in param_cov
param, param_cov = curve_fit(test, x, y)


print("Sine function coefficients:")
print(param)

print("Covariance of coefficients:")
print(param_cov)

np.sqrt(np.diag(param_cov))


# ans stores the new y-data according to
# the coefficients given by curve-fit() function
ans = (param[0]*(np.sin(param[1]*x)))

'''Below 4 lines can be un-commented for plotting results
using matplotlib as shown in the first example. '''

# plt.plot(x, y, 'o', color ='red', label ="data")
# plt.plot(x, ans, '--', color ='blue', label ="optimized data")
# plt.legend()
# plt.show()
