import numpy as np
import matplotlib.pyplot as plt


def mod(x):
    if x > 0 :
        return x
    else :
        return(-x)

def f_x(t) :
    return eta*Vlmax*t - (eta*m/Fmax)*(Vlmax**2)*(1-np.exp(-Ce*t))

def f_v(t) :
    return eta*Vlmax*(1-np.exp(-Ce*t))

def t0(s):
    a1 = eta*Vlmax
    b1 = (eta*m/Fmax)*(Vlmax**2)
    c1 = Ce
    
    t = (s/a1)
    for i in range(30):
        ft = a1*t - b1*(1-np.exp(-c1*t)) - s
        fdt = a1 - b1*c1*np.exp(-c1*t)
        h = -ft/fdt
        
        t = t + h
        
    return t
def r_best_cm(t):
    A = f/Tst
    G = (t*Tst)/(m*Wmax)
    
    rin = .1     # 5 cm
    for i in range(3000):
        Y = np.exp(-G/(rin**2))
        a11 = (Y-1)*2*A
        b11 = (Y-1)
        c11 = 2*(G**2)*(1-A*Y)
        
        d11 = (8*A*(G**2)*Y*(Y-1))
        e11 = 2*(G**2)*Y
        f11 = (4*A*(G**4)*Y)
        
        fr = a11*(rin**2) - b11*rin + c11
        fdr = (-d11/(rin**2)) + (1-Y) + (e11/(rin**2)) + (f11/(rin**3))
        
        h = -fdr/fr
        
        rin = rin + h
    return rin*100



def velocity_r(rin):
    return (1 - (rin*f/Tst))*(rin)*Wmax*(1 -np.exp((-G/(rin**2))))


def Rmax_Vmax():
    R = np.linspace(0.001,.25,2500)
    V = []
    for rin in R:
        vel = velocity_r(rin)
        V.append(vel)
    
    Vmxr = max(V)
    Rmax = 100*R[V.index(Vmxr)]     # cm
    rmx = Rmax/100
    return Rmax,Vmxr

gR = np.linspace(.24,2,80)
VelM = []
for gearR in gR:
   m = 60  # kg
   wmax_rpm = 6400*gearR
   Tst = 3.35/gearR  # Nm
   r_cm = 5    # cm
   mu = 0.02
   g = 10
   r = r_cm/100
   Wmax = wmax_rpm*(np.pi/30)
   rmax_cm = (Tst/(mu*m*g))*100
   Fmax = Tst/r
   f = mu*m*g
   Vlmax = r*Wmax
   Ce = Fmax/(m*Vlmax)
   eta = (1-(f/Fmax))
   ti = t0(850)
   A = f/Tst
   G = (ti*Tst)/(m*Wmax)
   
   VelocityM = Rmax_Vmax()[1]
   VelM.append(VelocityM)

plt.plot(gR,VelM)
plt.xlabel('Gear Ratio')
plt.ylabel('Velocity in m/s')
plt.show
    