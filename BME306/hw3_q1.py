import numpy as np
import matplotlib.pyplot as plt

# defining given constants
ker = 0.0284    # min^-1
kg = 0.0303     # min^-1
kpm = 0.00255   # min^-1
mu = 20

# cyto-b constants
ker2 = 0.0227    # min^-1
kg2 = 0.0192     # min^-1
kpm2 = 0.00612   # min^-1

# defining time bounds and empty arrays
t0 = 0
tEnd = 1000  # min
dt = 0.1
t = np.arange(t0, tEnd, 0.1)
n = int((tEnd - t0)/dt)
m = np.zeros(n)
per, pg, ppm, total = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
per2, pg2, ppm2, total2 = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)

# initializing arrays at t=0
per[0], total[0] = 20, 20
pg[0], ppm[0] = 0, 0

per2[0], total2[0] = 20, 20
pg2[0], ppm2[0] = 0, 0

# forward euler with differential equations 
for i in range(0, len(t)-1):
    per[i+1] = per[i] + dt * (-ker * per[i])
    pg[i+1] = pg[i] + dt * (-kg * pg[i] + ker*per[i])
    ppm[i+1] = ppm[i] + dt * (-kpm * ppm[i] + kg*pg[i])
    total [i+1] = per[i+1] + pg[i+1] + ppm[i+1]

    per2[i+1] = per2[i] + dt * (-ker2 * per2[i])
    pg2[i+1] = pg2[i] + dt * (-kg2 * pg2[i] + ker2*per2[i])
    ppm2[i+1] = ppm2[i] + dt * (-kpm2 * ppm2[i] + kg2*pg2[i])
    total2 [i+1] = per2[i+1] + pg2[i+1] + ppm2[i+1]



# plotting 
fig = plt.figure(num = 2, clear = True)
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, per, '-b', label = 'ER')
ax.plot(t, pg, '-r', label = 'Golgi')
ax.plot(t, ppm, '-g', label = 'PM')
ax.plot(t, total, label = 'Total')

ax.plot(t, per2, '--b', label = 'ER, Cyto-B')
ax.plot(t, pg2, '--r', label = 'Golgi, Cyto-B')
ax.plot(t, ppm2, '--g', label = 'PM, Cyto-B')
ax.plot(t, total2, linestyle = 'dashed', label = 'Total, Cyto-B')
ax.set(xlabel = 'time (minutes)', ylabel = 'VSVG-GfP (x10^6)', title = 'Cargo Protein (VSVG-GFP) in Compartments')
ax.grid(True), fig.tight_layout(), ax.legend()
plt.show()