from matplotlib import backend_bases
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression as LR
import numpy as np
from scipy.optimize import curve_fit

t = [0.5, 1.0, 2.0, 3.0, 5.0, 10, 20, 25, 30]
frac = [0.026, 0.052, 0.19, 0.23, 0.45, 0.69, 0.89, 0.81, 0.99]

def obj(t1, yr, tau):
    return 1-np.exp(-1 * (yr + (np.log(2)/15)) * (t1- tau))

params, pcov = curve_fit(obj, t, frac)

print(params)
#print("Optimal values for:\nGamma_r: ", params[0], "\nGrowth Factor, u: ", params[1], "\nDelay, tau: ", params[2])

fig = plt.figure(num = 2, clear = True)
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, obj(t, *params), '-r')
ax.scatter(t, frac)
ax.set(xlabel = 'time (hours)', ylabel = 'Fractional Radiolabeling (%)', title = 'Estimating half-life of polyA-mRNA')
ax.grid(True), fig.tight_layout()


plt.show()