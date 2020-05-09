import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def x(t):
    to = (Vmax)/(1-eta_a)
    t1 = t*(s-1)/2
    t2 = (1 + s + s*np.exp(s*t/tc) - np.exp(s*t/tc))/(2*s)
    
    return to*(t1 - tc*np.log(t2))

def v(t):
    ch = 1/np.arctanh(t/tca)
    den = 1 - s*ch
    
    return 2*Vmax*eta/den

def sol_S(S):
    ti = tc
    for t in range(50):
        h = -(x(ti) - S)/v(ti)
        ti = ti + h
    return ti


m = 60          # Infi Alpha Assumption
Tst =3.35       # Maxon BLDC
wrpm = 6400     # Maxon BLDC
Wmax = wrpm*np.pi/30
km = Tst/Wmax
Pmax = Tst*Wmax/4
g = 9.8
mu = 0.02       # assumption 

rcm = 6         # initial assumption which can be further optimized and studied
r = rcm/100

Vmax = r*Wmax
Fmax = Tst/r
f = mu*m*g
eta = (Fmax - f)/Fmax   # Friction Efficiency
tc = (m*Vmax)/Fmax      # Time constant Time to reach 63% of the maximum speed

Af = 0.1                # (1/2)*ro*Cd*A   = Fdrag/V^2
Fdmax = Af*(Vmax**2)    # Maximum Deag Force
eta_a = (Fmax - Fdmax)/Fmax     # Drag Efficiency

s = np.sqrt(1 + 4*eta*(1-eta_a))
B = Vmax/(1-eta_a)
D = B*s
tca = tc*(2/s)
T = np.linspace(0,150,1000)  # Time upto the pod reaches maximum acceleration distance which is 850 m
V = []
X = []
for t in T:
    velocity = v(t)
    pos = x(t)
    V.append(velocity)
    X.append(pos)

print(V[-1], "  The maximum velocity reached by the pod at end of 850 m")
print(T[-1], "seconds to reach the mark of 850m")
print(T[-1]/tc, " ratio of time taken to reach 850m to time constant")
plt.plot(T,V)
plt.xlabel("Time in s")
plt.ylabel("Velocity in m/s")
plt.grid(True)

 