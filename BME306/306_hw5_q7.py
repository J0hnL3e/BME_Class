import numpy as np
import matplotlib.pyplot as plt

# defining given constants
k1 = 0.143  # day^-1
k2 = 1      # day^-1
k3 = 0.01   # day^-1
k4 = 5.35   # day^-1
k5 = 30     # day^-1
k6 = 0.01   # day^-1
k7 = 0.1    # day^-1
k8 = 1      # day^-1

# modified k1 rates 
k1s = 0.2 * k1
k1m = 0.15 * k1
k1l = 0.1 * k1

t0 = 0
tEnd = 200  # day
dt = 0.01
t = np.arange(t0, tEnd, 0.01)
n = int((tEnd - t0)/dt)

p, p2, p3, p4 = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
csc, csc2, csc3, csc4 = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
d, d2, d3, d4 = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
total, total2, total3, total4 = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
N = 1e4  # cells
vc = 4.18*1e-6
p[0], d[0], csc[0] = (0.1991 * N), (0.7991 * N), (0.0018 * N)
p2[0], d2[0], csc2[0] = (0.1991 * N), (0.7991 * N), (0.0018 * N)
p3[0], d3[0], csc3[0] = (0.1991 * N), (0.7991 * N), (0.0018 * N)
p4[0], d4[0], csc4[0] = (0.1991 * N), (0.7991 * N), (0.0018 * N)
total[0], total2[0], total3[0], total4[0] = p[0] + d[0] + csc[0], p[0] + d[0] + csc[0], p[0] + d[0] + csc[0], p[0] + d[0] + csc[0]

for i in range(0, len(t)-1):
    if (t[i] == 120):
        p[i+1] = p[i] + dt * (csc[i] * (k2 + 2*k3) + p[i] * (k4-k5-k8))
        d[i+1] = d[i] + dt * (2*k5*p[i]-k7*d[i])
        csc[i+1] = csc[i] + dt * (csc[i]*(k1-k3-k6))
        total[i+1] = p[i+1] + d[i+1] + csc[i+1]

        p2[i+1] = 0.1 * p2[i]
        d2[i+1] = 0.1 * d2[i]
        csc2[i+1] = 0.1 * csc2[i]
        total2[i+1] = p2[i+1] + d2[i+1] + csc2[i+1]

        p3[i+1] = 0.01 * p3[i]
        d3[i+1] = 0.01 * d3[i]
        csc3[i+1] = 0.01 * csc3[i]
        total3[i+1] = p3[i+1] + d3[i+1] + csc3[i+1]

        p4[i+1] = 0.001 * p4[i]
        d4[i+1] = 0.001 * d4[i]
        csc4[i+1] = 0.001 * csc4[i]
        total4[i+1] = p4[i+1] + d4[i+1] + csc4[i+1]

    else:
        p[i+1] = p[i] + dt * (csc[i] * (k2 + 2*k3) + p[i] * (k4-k5-k8))
        d[i+1] = d[i] + dt * (2*k5*p[i]-k7*d[i])
        csc[i+1] = csc[i] + dt * (csc[i]*(k1-k3-k6))
        total[i+1] = p[i+1] + d[i+1] + csc[i+1]

        p2[i+1] = p2[i] + dt * (csc2[i] * (k2 + 2*k3) + p2[i] * (k4-k5-k8))
        d2[i+1] = d2[i] + dt * (2*k5*p2[i]-k7*d2[i])
        csc2[i+1] = csc2[i] + dt * (csc2[i]*(k1-k3-k6))
        total2[i+1] = p2[i+1] + d2[i+1] + csc2[i+1]

        p3[i+1] = p3[i] + dt * (csc3[i] * (k2 + 2*k3) + p3[i] * (k4-k5-k8))
        d3[i+1] = d3[i] + dt * (2*k5*p3[i]-k7*d3[i])
        csc3[i+1] = csc3[i] + dt * (csc3[i]*(k1-k3-k6))
        total3[i+1] = p3[i+1] + d3[i+1] + csc3[i+1]

        p4[i+1] = p4[i] + dt * (csc4[i] * (k2 + 2*k3) + p4[i] * (k4-k5-k8))
        d4[i+1] = d4[i] + dt * (2*k5*p4[i]-k7*d4[i])
        csc4[i+1] = csc4[i] + dt * (csc4[i]*(k1-k3-k6))
        total4[i+1] = p4[i+1] + d4[i+1] + csc4[i+1]

res = [vc * i for i in total]
res2 = [vc * i for i in total2]
res3 = [vc * i for i in total3]
res4 = [vc * i for i in total4]

# plotting 
fig = plt.figure(num = 1, clear = True)
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, res, label = 'Tumor')
ax.set(xlabel = 'time (days)', ylabel = 'Tumor Volume (mm^3)')
ax.grid(True), fig.tight_layout(), ax.legend()
plt.xlim([50, 150]), plt.ylim([0, 1400])

fig = plt.figure(num = 2, clear = True)
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, res, label = 'Tumor, no non-target chemo')
ax.plot(t, res2, label = 'Tumor, 90% non-target chemo')
ax.plot(t, res3, label = 'Tumor, 99% non-target chemo')
ax.plot(t, res4, label = 'Tumor, 99.9% non-target chemo')
ax.set(xlabel = 'time (days)', ylabel = 'Tumor Volume (mm^3)')
ax.grid(True), fig.tight_layout(), ax.legend()
plt.ylim([0, 6000])
plt.show()