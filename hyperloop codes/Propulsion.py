import numpy as np
import matplotlib.pyplot as plt


m = 60  # kg
gearR = 2
wmax_rpm = 6400*gearR
Tst = 3.35/gearR  # Nm
r_cm = 3.75    # cm
mu = 0.02
g = 10

r = r_cm/100
Wmax = wmax_rpm*(np.pi/30)


Fmax = Tst/r
f = mu*m*g
Vlmax = r*Wmax
Ce = Fmax/(m*Vlmax)

eta = (1-(f/Fmax))


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



ti = t0(850)
A = f/Tst
G = (ti*Tst)/(m*Wmax)

R = np.linspace(0.01,.2,1000)
V = []
for rin in R:
    vel = (1 - (rin*f/Tst))*(rin)*Wmax*(1 -np.exp((-G/(rin**2))))
    V.append(vel)
    
Vmxr = max(V)
Rmax = 100*R[V.index(Vmxr)]     # cm
rmx = Rmax/100
    
plt.plot(100*R,V)
plt.xlabel('Radius in cm')
plt.ylabel('Velocity in m/s')
plt.grid()
plt.show    


#print(f_v(t0(850)),r_best_cm(t0(850)))

#gr = np.linspace(1/3,3,1000)
#Vmaxgr = []
#for gearR in gr:
#    Tst = Tst/gearR
#    Wmax = Wmax*gearR
#    
#    vel = (1 - (rmx*f/Tst))*(rmx)*Wmax*(1 -np.exp((-G/(rmx**2))))
#    Vmaxgr.append(vel)
#
#plt.plot(gr,Vmaxgr)   
    