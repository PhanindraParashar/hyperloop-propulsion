import numpy as np
import matplotlib.pyplot as plt

m = .6          # Infi Alpha Assumption
n = 1

tst_o = .41
w_m_o = 16800*(np.pi/30)

Tst = tst_o/n     # Maxon BLDC

Wmax = n*w_m_o
kmo = tst_o/w_m_o

km = kmo/(n**2)
Pmax = Tst*Wmax/4
g = 9.8
mu = 0.05       # assumption 
f = mu*m*g

#rcm =  100*(tst_o/(2*f))*1        # initial assumption which can be further optimized and studied
rcm = 3.5
r = rcm/100

Vmax = r*Wmax
Fmax = Tst/r

eta = (Fmax - f)/Fmax   # Friction Efficiency
tc = (m*Vmax)/Fmax      # Time constant Time to reach 63% of the maximum speed

Af = 0.0005            # (1/2)*ro*Cd*A   = Fdrag/V^2
Fdmax = Af*(Vmax**2)    # Maximum Deag Force
eta_a = (Fmax - Fdmax)/Fmax     # Drag Efficiency

s = np.sqrt(1 + 4*eta*(1-eta_a))
B = Vmax/(1-eta_a)
C = B*eta*Vmax


def v(t):
    to = eta*Vmax
    t1 = 1 - np.exp(-s*(t/tc))
    t2 = s + 1 + (s-1)*np.exp(-s*(t/tc))
    return to*(t1/t2)



def x(t):
    to = Vmax/(4*(1-eta_a))
    t1 = t*(s-1)
    t3 = (s + 1 + (s-1)*np.exp(-s*(t/tc)))/(2*s)
    return to*(t1 + 2*tc*np.log(t3))

def sol_x(s):   # solving for the equation to get the time when pod reaches distance s.  Newtons approximation method is applied 
    ti = tc     # initial approximation, it could be anything. taking t = tc for simplicity
    for i in range(50): # while loop with condition on accurecy can also be used, however even 5 to 6 iterations generally give very good approximation. In this case we are doing 50 iterations!!!
        h = -(x(ti) - s)/v(ti)
        ti = ti + h
    return ti 



T = np.linspace(0,sol_x(100),1000)  # Time upto the pod reaches maximum acceleration distance which is 850 m
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
print(X[-1], " distance in m")
print(T[-1]/tc, " ratio of time taken to reach ___m to time constant")
plt.plot(T,V)
plt.xlabel("Time in s")
plt.ylabel("Velocity in m/s")
plt.grid(True)