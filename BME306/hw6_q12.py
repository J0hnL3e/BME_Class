import numpy as np
import matplotlib.pyplot as plt

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

gluc, amm, lac, igG = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
v, d, n, glut = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)

for i in range(0, len(t)-1): 
    rp = rpm * (glut[i]/(Kglu + glut[i])) * (k_ai/(k_ai+amm[i])) * (k_li/(k_li+lac[i]))
    rd = rdm * (amm[i]/(Ka + amm[i])) * (lac[i]/(Kl+lac[i])) * (K_glui/(K_glui+glut[i]))

    v[i+1] = v[i] + dt * (rp - rd) * v[i]
    d[i+1] = d[i] + dt * (rd * v[i])
    n[i+1] = v[i+1] + d[i+1]
    gluc[i+1] = gluc[i] + dt * ((-1/yg) * rp * v[i] - m*v[i])

    amm[i+1] = (-1/ya) * 
    dglut = (-1/yglu) * rp * v[i] - (1/yp) * () - (1/yag) * ()
    glut[i+1] = 
    lac[i+1] = lac[i] + dt * (-1/yl) * ((-1/yg) * rp * v[i] - m*v[i])



