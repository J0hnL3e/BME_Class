import numpy as np
import matplotlib.pyplot as plt

ker = 0.0284    # min^-1
kg = 0.0303     # min^-1
kpm = 0.00255   # min^-1
mu = 20

t0 = 0
tEnd = 600  # min
dt = 0.1

t = np.arange(t0, tEnd, 0.1)
n = int((tEnd - t0)/dt)
m = np.zeros(n)
per = np.zeros(n)
pg = np.zeros(n)
ppm = np.zeros(n)

per[0] = 20
pg[0], ppm[0] = 0, 0

for i in range(0, len(t)-1):
    per[i+1] = per[i] + dt * (-ker * per[i])
    pg[i+1] = pg[i] + dt * (-kg * pg[i] + ker*per[i])
    ppm[i+1] = ppm[i] + dt * (-kpm * ppm[i] + kg*pg[i])

fig = plt.figure(num = 2, clear = True)
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, per, '-b', label = 'Analyltic mRNA')
ax.plot(t, pg, '-r', label = 'mRNA_ss constants')
ax.plot(t, ppm, '--g', label = 'dt = 0.1')
ax.set(xlabel = 'time (hours)', ylabel = 'Proteins (molecule/uM^3)', title = 'Comparison of Analytical and Numerical')
ax.grid(True), fig.tight_layout(), ax.legend()
plt.show()