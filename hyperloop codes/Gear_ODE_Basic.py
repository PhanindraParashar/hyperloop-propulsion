import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

n = .7
m = 60          # Infi Alpha Assumption
def get_constants(n):
    Tst = 24.12/n       # Maxon BLDC
    wrpm = 4400*n     # Maxon BLDC
    Wmax = wrpm*np.pi/30
    km = Tst/Wmax
    Pmax = Tst*Wmax/4
    g = 9.8
    mu = 0.01       # assumption 
    
    rcm = 14         # initial assumption which can be further optimized and studied
    r = rcm/100
    
    Vmax = r*Wmax
    Fmax = Tst/r
    f = mu*m*g
    eta = (Fmax - f)/Fmax   # Friction Efficiency
    tc = (m*Vmax)/Fmax      # Time constant Time to reach 63% of the maximum speed
    
    Af = 0.003                # (1/2)*ro*Cd*A   = Fdrag/V^2
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
    
    return C1,C2,C3


def Graphs(n,tmx,x0,v0):
    # initial condition
    y0 = [x0,v0]
    
    # time points
    t = np.linspace(0,tmx,1000)
    
    # solve ODE
    sol = odeint(derivatives, y0, t)
    x = sol[:, 0]
    v = sol[:, 1]
    
    xl = x[-1]
    vl = v[-1]
    # plot results
    
    #plt.plot(t, sol[:, 0], 'b', label='x(t)')
#    plt.plot(t, sol[:, 1], 'g', label='v(t)')
#    
#    plt.xlabel('t')
#    plt.grid()
#    plt.show()
    
    return xl,vl
    

def model(v,t): 
    (C1,C2,C3) = get_constants(n)
    dvdt = C1 - C2*v - C3*(v**2)
    return dvdt

def derivatives(y, t):
    (C1,C2,C3) = get_constants(n)
    x, v = y
    dydt = [v, C1 - C2*v - C3*(v**2)]
    return dydt

nst = .5
d = .5
nt = 3
ni = np.linspace(nst,nst + nt*d,nt+1)
#ni = [nst,nst+d,nst+2*d,nst+3*d,nst+4*d,nst + 5*d,nst +6*d]

tstep = 32/len(ni)
x0 = 0
v0 = 0
for n in ni:
    
    (x0,v0) = Graphs(n,tstep,x0,v0)
    print(x0,v0)
    tstep = tstep  




