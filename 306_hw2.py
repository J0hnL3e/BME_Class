#from sklearn.linear_model import LinearRegression as LR
import numpy as np
import matplotlib.pyplot as plt

# Given Values
kr = 0.0015 # molecules of mRNA /(h * um)^3
yr = 0.075  # hr^-1
kp = 10     # molecules of Protein/(molecules of mRNA*hour)
yp = 0.01   # hr^-1
u = 0.05    # hr^-1
t0 = 0
tEnd = 100  # hours
# dt = 0.1

def analytic():
    # using mRNA, analytic
    anal1 = [((kp*kr)/((yp + u) * (yr + u)))\
        *(1-np.exp((-yr -u)*i))\
        *(1+(  (yr+u)*(np.exp(-1*(yp+u)*i)) - (yp+u)*(np.exp(-1*(yr+u)*i))  )/(yp-yr)) for i in t]

    # just using mRNA_ss
    anal_ss = [((kp*kr)/((yp + u) * (yr + u)))\
        *(1+(  (yr+u)*(np.exp(-1*(yp+u)*i)) - (yp+u)*(np.exp(-1*(yr+u)*i))  )/(yp-yr)) for i in t]
    
    return anal1, anal_ss

t = np.arange(t0, tEnd, 0.1)
# t5 = np.arange(t0, tEnd, 5)
# t10 = np.arange(t0, tEnd, 10)
# n = int((tEnd - t0)/dt)
# m = np.zeros(n)
# p1 = np.zeros(n)
# p2 = np.zeros(int((tEnd - t0)/5))
# p3 = np.zeros(int((tEnd - t0)/10))
# print(len(t), len(t5), len(t10))
# print(len(p1), len(p2), len(p3))

# m[0]  = 0
# p1[0] = 0
# p2[0] = 0
# p3[0] = 0

# for i in range (0, len(t)-1):
#     m[i+1] = m[i] + dt * (kr - m[i]*(yr + u))
#     p1[i+1] = p1[i] + dt * (kp * m[i] - p1[i] * (yp + u))
# for j in range(0, len(t5)-1):
#     p2[i+1] = p2[i] + 5 * (kp * m[i] - p2[i] * (yp + u))
# for k in range(0, len(t10)-1):
#     p3[i+1] = p3[i] + 10 * (kp * m[i] - p3[i] * (yp + u))

#x = [item*500*np.exp(2.3) for item in t]

def num(dt):
    t0, tEnd = 0, 100
    n = int((tEnd-t0)/dt)
    t = np.arange(0, 100, dt)
    m = np.zeros(n)
    p = np.zeros(n)
    m[0], p[0] = 0, 0

    for i in range (0, len(t)-1):
        m[i+1] = m[i] + dt * (kr - m[i]*(yr + u))
        p[i+1] = p[i] + dt * (kp * m[i] - p[i] * (yp + u))
    
    return t, m, p

fig = plt.figure(num=1, clear=True)
ax = fig.add_subplot(1,1,1)
ax.plot(t, num(0.1)[2], '-g', label = "dt = 0.1")
ax.plot(num(5)[0], num(5)[2], '-b', label = "dt = 5")
ax.plot(num(10)[0], num(10)[2], '-r', label = "dt = 10")
# ax.plot(t5, p2, label = "dt = 5")
# ax.plot(t10, p3, label = "dt = 10")
ax.set(xlabel = 'time (hours)', ylabel = 'Proteins (molecule count)', title = 'Comparison of Analytical and Numerical')
ax.grid(True), fig.tight_layout(), ax.legend()

# Testing out Different [mRNA] for [P]ss
fig = plt.figure(num = 2, clear = True)
ax = fig.add_subplot(1, 1, 1)
ax.plot(t, analytic()[0], '-b', label = 'Analyltic mRNA')
ax.plot(t, analytic()[1], '-r', label = 'mRNA_ss constants')
ax.plot(t, num(0.1)[2], '--g', label = 'dt = 0.1')
ax.set(xlabel = 'time (hours)', ylabel = 'Proteins (molecule count)', title = 'Comparison of Analytical and Numerical')
ax.grid(True), fig.tight_layout(), ax.legend()

plt.show()

