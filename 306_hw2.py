from sklearn.linear_model import LinearRegression as LR
import numpy as np
import math
import matplotlib.pyplot as plt

t0 = 0
tEnd = 100  # hours
dt = 0.1
t = np.arange(t0, tEnd, dt)
t5 = np.arange(t0, tEnd, 5)
t10 = np.arange(t0, tEnd, 10)
n = int((tEnd - t0)/dt)

m = np.zeros(n)
p1 = np.zeros(n)
p2 = np.zeros(n)
p3 = np.zeros(n)

m[0]  = 0
p1[0] = 0
p2[0] = 0
p3[0] = 0

kr = 0.0015 # molecules of mRNA /(h * um)^3
yr = 0.075  # hr^-1
kp = 10     # molecules of Protein/(molecules of mRNA*hour)

yp = 0.01   # hr^-1
u = 0.05    # hr^-1

# using mRNA, analytic
protein_analytic = [((kp*kr)/((yp + u) * (yr + u)))\
    *(1-np.exp((-yr -u)*i))\
    *(1+(  (yr+u)*(np.exp(-1*(yp+u)*i)) - (yp+u)*(np.exp(-1*(yr+u)*i))  )/(yp-yr)) for i in t]

# just using mRNA_ss
protein_analytic_2 = [((kp*kr)/((yp + u) * (yr + u)))\
    *(1+(  (yr+u)*(np.exp(-1*(yp+u)*i)) - (yp+u)*(np.exp(-1*(yr+u)*i))  )/(yp-yr)) for i in t]

for i in range (0, len(t)-1):
    m[i+1] = m[i] + dt * (kr - m[i]*(yr + u))
    p1[i+1] = p1[i] + dt * (kp * m[i] - p1[i] * (yp + u))

    if i % 50 == 0:
        p2[i+1] = p2[i] + 5 * (kp * m[i] - p2[i] * (yp + u))
    
    else
    

    
    if i % 10 == 0:
        p3[i+1] = p3[i] + 10 * (kp * m[i] - p3[i] * (yp + u))


#x = [item*500*np.exp(2.3) for item in t]
x = [item for item in t]

fig = plt.figure(num=1, clear=True)
ax = fig.add_subplot(1,1,1)
ax.plot(t, protein_analytic, '-b', label = "Analytic mRNA")
#ax.plot(t, protein_analytic_2, '-r', label = "mRNA_ss")
# ax.plot(t, p1, '--g')
ax.plot(t, p1, '-g', label = "dt = 0.1")
ax.plot(t, p2, label = "dt = 5")
ax.plot(t, p3, label = "dt = 10")
ax.set(xlabel = 'time (hours)', ylabel = 'Proteins (molecule count)', title = 'Comparison of Analytical and Numerical')
ax.grid(True), fig.tight_layout(), ax.legend()

# fig = plt.figure(num=2, clear=True)
# ax = fig.add_subplot(1,1,1)
# ax.plot(t, m, '-b')
# ax.plot(t, p1, '-r')
# ax.set(xlabel = 'time (hours)', ylabel = 'Proteins (molecule count)', title = 'Comparison of Analytical and Numerical')
# ax.grid(True), fig.tight_layout()
plt.show()

print(n)
print(len(t5))