import numpy as np
import matplotlib.pyplot as plt

# defining given constants
ker = 0.0284    # min^-1
kg = 0.0303     # min^-1
kpm = 0.00255   # min^-1
mu = 20

# defining time bounds and empty arrays
t0 = 0
tEnd = 600  # min
dt = 0.1
t = np.arange(t0, tEnd, 0.1)
n = int((tEnd - t0)/dt)
m = np.zeros(n)
per = np.zeros(n)
pg = np.zeros(n)
ppm = np.zeros(n)
total= np.zeros(n)

# initializing arrays at t=0
per[0] = 20
pg[0], ppm[0], total[0] = 0, 0, 20

# forward euler with differential equations
for i in range(0, len(t)-1):
    per[i+1] = per[i] + dt * (-ker * per[i])
    pg[i+1] = pg[i] + dt * (-kg * pg[i] + ker*per[i])
    ppm[i+1] = ppm[i] + dt * (-kpm * ppm[i] + kg*pg[i])
    total [i+1] = per[i+1] + pg[i+1] + ppm[i+1]

# plotting 
fig = plt.figure(num = 2, clear = True)
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, per, '-b', label = 'ER')
ax.plot(t, pg, '-r', label = 'PM')
ax.plot(t, ppm, '-g', label = 'Golgi')
ax.plot(t, total, label = 'Total')
ax.set(xlabel = 'time (minutes)', ylabel = 'VSVG-GfP (x10^6)', title = 'Cargo Protein (VSVG-GFP) in Compartments')
ax.grid(True), fig.tight_layout(), ax.legend()
plt.show()