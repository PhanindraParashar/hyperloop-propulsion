import numpy as np
import matplotlib.pyplot as plt

m = 60          # Infi Alpha Assumption
n = 1
Tst = 24.12/n       # Maxon BLDC
wrpm = 4400*n     # Maxon BLDC
Wmax = wrpm*np.pi/30
km = Tst/Wmax
Pmax = Tst*Wmax/4
g = 9.8
mu = 0.05       # assumption 
f = mu*m*g

rcm =  100*(Tst/(2*f))        # initial assumption which can be further optimized and studied
#rcm = 20
r = rcm/100

Vmax = r*Wmax
Fmax = Tst/r

eta = (Fmax - f)/Fmax   # Friction Efficiency
tc = (m*Vmax)/Fmax      # Time constant Time to reach 63% of the maximum speed

def x(t):       # function that returns distance on input of time,  this is derived from solving differential equation under initial condition at t = 0   x = 0
    return eta*Vmax*(t -tc + tc*np.exp(-t/tc))

def v(t):       # function that is derived from solving differential equation and condition that at t = 0   v = 0
    return eta*Vmax*(1-np.exp(-t/tc))

def a(t):       # function that returns acceleration
    return (eta*Vmax/tc)*np.exp(-t/tc)

def sol_x(s):   # solving for the equation to get the time when pod reaches distance s.  Newtons approximation method is applied 
    ti = tc     # initial approximation, it could be anything. taking t = tc for simplicity
    for i in range(50): # while loop with condition on accurecy can also be used, however even 5 to 6 iterations generally give very good approximation. In this case we are doing 50 iterations!!!
        h = -(x(ti) - s)/v(ti)
        ti = ti + h
    return ti 

# Plotting velocity position and Time
T = np.linspace(0,4*tc,1000)  # Time upto the pod reaches maximum acceleration distance which is 850 m
V = []
X = []
for t in T:
    velocity = v(t)
    pos = x(t)
    V.append(velocity)
    X.append(pos)

print(rcm, " radius in cm")
print(V[-1], "  The maximum velocity reached by the pod at end of graph")
print(T[-1], "seconds to reach")
print(T[-1]/tc, " ratio of time taken to reach ___m to time constant")
plt.plot(T,V)
plt.xlabel("Time in s")
plt.ylabel("Velocity in m/s")
plt.grid(True)