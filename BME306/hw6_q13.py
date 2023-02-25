from matplotlib import backend_bases
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression as LR
import numpy as np
from scipy.optimize import curve_fit

t0 = 0
tEnd = 20  # concentration (g/L)
dt = 0.01
t = np.arange(t0, tEnd, dt)
tp = np.arange(0, 5, 1)
n = int((tEnd - t0)/dt)
gluc = [0.667, 1.124, 2.353, 3.617, 4.886]
prof = [0.005, 0.019, 0.028, 0.033, 0.036]

def obj(x, rm, ks):
    y = rm * (x/(ks+x))
    return y

params, pcov = curve_fit(obj, gluc, prof)
print(params)
#print("Optimal values for:\nGamma_r: ", params[0], "\nGrowth Factor, u: ", params[1], "\nDelay, tau: ", params[2])
ci = 2.776 * np.std(prof)/np.sqrt(len(gluc))
# df = 5-1 = 4
# due to small sample size, t-value is 2.776 for 95% 

new_fit = obj(gluc, *params)
ints = []
for i in params:
    f = (i-ci, i+ci)
    ints.append(f)
print(ints)

minConc = obj(t, *params)
for i in range(0, len(minConc) - 1):
    if (minConc[i] > 0):
        mmm = t[i]
        break
print(mmm)

fig = plt.figure(num = 3, clear = True)
ax = fig.add_subplot(1, 1, 1)
ax.plot(gluc, new_fit, '-r')
ax.scatter(gluc, prof)
ax.fill_between(gluc, (new_fit-ci), (new_fit+ci), color='b', alpha=.1)
ax.set(xlabel = 'Glucose Conc. (g/L)', ylabel = 'Proliferation Rate (h^-1)', title = 'Monod Kinetic Proliferation Model')
ax.grid(True), fig.tight_layout()

fig = plt.figure(num = 4, clear = True)
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, obj(t, *params), '-r')
ax.scatter(gluc, prof)
ax.fill_between(gluc, (new_fit-ci), (new_fit+ci), color='b', alpha=.1)
ax.set(xlabel = 'Glucose Conc. (g/L)', ylabel = 'Proliferation Rate (h^-1)', title = 'Monod Kinetic Proliferation Model, Extended')
ax.grid(True), fig.tight_layout()
plt.show()

