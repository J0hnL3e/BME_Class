import numpy as np
import matplotlib.pyplot as plt
import itertools

rpm, rdm, ksec =  3, 2, 0.8 # day^-1
riggm = 1.5e-7  # mg IgG/(mL*day*cells/mL)
m = 1.25e-8     # mM
Kglu, Ka, Kl = 0.8, 1.44, 15 # mM
K_glui, k_ai, k_li, k_ni = 5e-4, 1.05, 8, 1e7 # mM
kglup = 0.01 # mM
yg, yglu = 2.9e5, 9.434e5 
yp, yag, ya, yl = 1000, 2.1937, 1.5, 3.00 # mM/mM

t0 = 0
tEnd = 30  # days
dt = 0.01
t = np.arange(t0, tEnd, dt)
n = int((tEnd - t0)/dt)

gluc, amm, lac, igG, egg = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
v, d, n, glut = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)

n[0] = 30e4
d[0] = 0.05 * n[0]
v[0] = 0.95 * n[0]
gluc[0] = 22.3
glut[0] = 2.5
amm[0] = 0.8
lac[0] = 3.25
egg[0] = 0

for i in range(0, len(t)-1): 
    rp = rpm * (glut[i]/(Kglu + glut[i])) * (k_ai/(k_ai+amm[i])) * (k_li/(k_li+lac[i]))
    rd = rdm * (amm[i]/(Ka + amm[i])) * (lac[i]/(Kl+lac[i])) * (K_glui/(K_glui+glut[i]))

    v[i+1] = v[i] + dt * (rp - rd) * v[i]
    d[i+1] = d[i] + dt * (rd * v[i])
    n[i+1] = v[i+1] + d[i+1]
    gluc[i+1] = gluc[i] + dt * ((-1/yg) * rp * v[i] - m*v[i])
    lac[i+1] = lac[i] + dt * (-1/yl) * ((-1/yg) * rp * v[i] - m*v[i])
    digg = riggm * (glut[i]/(kglup+glut[i])) * (n[i]/(k_ni+n[i])) * v[i]
    dglut = ((-1/yglu) * rp * v[i] - (1/yp) * (digg))/(1-1/(yag*ya))
    amm[i+1] = (-1/ya) * dglut * dt + amm[i]
    
    glut[i+1] = glut[i] + dglut * dt
    
    igG[i+1] = dt * digg + igG[i]
    egg[i+1] = dt * ksec *(igG[i]-egg[i]) + egg[i]

percs = []
for (i, j) in zip(v, n):
    percs.append(i/j)

#percs = [i/j for i in v and j in t]

fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, v, label = 'Viable Cells')
ax.plot(t, n, label = 'Total Cells')
ax.set(xlabel = 'time (days)', ylabel = 'Density (cells/mL)')
ax.grid(True), fig.tight_layout(), ax.legend()

fig = plt.figure(num = 2, clear = True)
ax1 = fig.add_subplot(3, 2, 1)
ax2 = fig.add_subplot(3, 2, 2)
ax3 = fig.add_subplot(3, 2, 3)
ax4 = fig.add_subplot(3, 2, 4)
ax5 = fig.add_subplot(3, 2, 5)
ax6 = fig.add_subplot(3, 2, 6)

ax1.plot(t, percs)
ax1.set(xlabel = 'time (days)', ylabel = 'Percent Viable (%)', title = 'Percent Cell Viability')
ax1.grid(True), fig.tight_layout()

ax2.plot(t, gluc)
ax2.set(xlabel = 'time (days)', ylabel = 'Glucose Conc. (mM)', title = 'Glucose')
ax2.grid(True), fig.tight_layout()

ax3.plot(t, glut)
ax3.set(xlabel = 'time (days)', ylabel = 'Glutamine Conc. (mM)', title = 'Glutamine')
ax3.grid(True), fig.tight_layout()

ax4.plot(t, amm)
ax4.set(xlabel = 'time (days)', ylabel = 'Ammonia Conc. (mM)', title = 'Ammonia')
ax4.grid(True), fig.tight_layout()

ax5.plot(t, lac)
ax5.set(xlabel = 'time (days)', ylabel = 'Lactate Conc. (mM)', title = 'Lactate')
ax5.grid(True), fig.tight_layout()

ax6.plot(t, egg)
ax6.set(xlabel = 'time (days)', ylabel = 'Extracellular igG Conc. (mg/L)', title = 'Extracellular igG')
ax6.grid(True), fig.tight_layout()
plt.show()




