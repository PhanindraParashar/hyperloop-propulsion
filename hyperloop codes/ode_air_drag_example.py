import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def model(v,t):    
    dvdt = C1 - C2*v - C3*(v**2)
    return dvdt

def derivatives(y, t):
    x, v = y
    dydt = [v, C1 - C2*v - C3*(v**2)]
    return dydt


def coth(x):
    return 1/np.tanh(x)



def R_max():
    R = np.linspace(0.01,2,100)
    Fmax = Tst/R
    E = (Fmax - f)/Fmax
    A = (Fmax - Fdmax)/Fmax
    Ve = []
    s = np.sqrt(1 + 4*E*A)
    for r in R:
        v = 2*E*r*Wmax/(1 - s)
        Ve.append(v)
    
    #vmx = max(Ve)
    #rind = Ve.any(index(vmx))
    #rmax = R[rind]
    return Ve#rmax


n = 1
m = 60          # Infi Alpha Assumption
Tst = 3.35       # Maxon BLDC
wrpm = 6400     # Maxon BLDC
Wmax = wrpm*np.pi/30
km = Tst/Wmax
Pmax = Tst*Wmax/4
g = 9.8
mu = 0.01       # assumption 
f = mu*m*g

rcm = 100*(Tst/(2*f))         # initial assumption which can be further optimized and studied

r = rcm/100

Vmax = r*Wmax
Fmax = Tst/r

eta = (Fmax - f)/Fmax   # Friction Efficiency
tc = (m*Vmax)/Fmax      # Time constant Time to reach 63% of the maximum speed

Af = 0.00005                # (1/2)*ro*Cd*A   = Fdrag/V^2
Fdmax = Af*(Vmax**2)    # Maximum Deag Force
eta_a = (Fmax - Fdmax)/Fmax     # Drag Efficiency

s = np.sqrt(1 + 4*eta*(1-eta_a))
B = Vmax/(1-eta_a)
C = B*eta*Vmax
D = B*s
tca = tc*(2/s)

C1 = (Tst/(m*r)) - mu*g
C2 = km/(m*r*r)
C3 = Af/m
    
# initial condition
y0 = [0,0]

# time points
t = np.linspace(0,5*tc,10000)

# solve ODE
sol = odeint(derivatives, y0, t)
x = sol[:, 0]
v = sol[:, 1]

print(x[-1])
print(v[-1])
print(rcm)
# plot results

#plt.plot(t, sol[:, 0], 'b', label='x(t)')
plt.plot(t, sol[:, 1], 'g', label='v(t)')

plt.xlabel('t')
plt.grid()
plt.show()
    
    










