import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

m = 130         # Infi Alpha Assumption
Tst = 1.65*2*2       # Maxon BLDC
wrpm = 6900/2     # Maxon BLDC
Wmax = wrpm*np.pi/30
km = Tst/Wmax
Pmax = Tst*Wmax/4
g = 9.8
mu = 0.02       # assumption 

rcm = 10         # initial assumption which can be further optimized and studied
r = rcm/100

Vmax = r*Wmax
Fmax = Tst/r
f = mu*m*g
eta = (Fmax - f)/Fmax   # Friction Efficiency
tc = (m*Vmax)/Fmax      # Time constant Time to reach 63% of the maximum speed

#Af = 0.1/20                # (1/2)*ro*Cd*A   = Fdrag/V^2
Af = 0.05
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

def model(v,t):    
    dvdt = C1 - C2*v - C3*(v**2)
    return dvdt

def derivatives(y, t, C1,C2,C3):
    x, v = y
    dydt = [v, C1 - C2*v - C3*(v**2)]
    return dydt




# initial condition
y0 = [0,0]

# time points
#tstop = 17.4
tstop = 3*tc
t = np.linspace(0,tstop,1000)

# solve ODE
sol = odeint(derivatives, y0, t, args=(C1,C2,C3))
x = sol[:, 0]/1000
v = sol[:, 1]

print(x[-1], " Km")
print(v[-1]*3.6, " Kmph")
print(t[-1]/60, " minutes")
# plot results

#plt.plot(t, sol[:, 0], 'b', label='x(t)')
plt.plot(t, v, 'g', label='v(t)')

plt.xlabel('t')
plt.grid()
plt.show()





