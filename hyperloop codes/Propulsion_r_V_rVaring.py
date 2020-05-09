import numpy as np
import matplotlib.pyplot as plt


def f_x(t,r) :
    Fmax = Tst/r
    f = mu*m*g
    Vlmax = r*Wmax
    Ce = Fmax/(m*Vlmax)
    
    eta = (1-(f/Fmax))
    return eta*Vlmax*t - (eta*m/Fmax)*(Vlmax**2)*(1-np.exp(-Ce*t))

def f_v(t,r) :
    Fmax = Tst/r
    f = mu*m*g
    Vlmax = r*Wmax
    Ce = Fmax/(m*Vlmax)
    
    eta = (1-(f/Fmax))
    return eta*Vlmax*(1-np.exp(-Ce*t))


def t0(s,r):
    Fmax = Tst/r
    f = mu*m*g
    Vlmax = r*Wmax
    Ce = Fmax/(m*Vlmax)

    eta = (1-(f/Fmax))
    
    
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


m = 60  # kg
gearR = 1
wmax_rpm = 6400*gearR
Tst = 3.35/gearR  # Nm

mu = 0.0
g = 10


Wmax = wmax_rpm*(np.pi/30)

km = Tst/Wmax


R = np.linspace(0.01,.12,100)
V = []
for rin in R:
    f = mu*m*g
    A = f/Tst
    ti = t0(850,rin)
    G = (ti*Tst)/(m*Wmax)
    
    vel = (1 - (rin*f/Tst))*(rin)*Wmax*(1 -np.exp((-G/(rin**2))))
    V.append(vel)
    
Vmxr = max(V)
Rmax = 100*R[V.index(Vmxr)]     # cm
rmx = Rmax/100

tc = m*rmx*rmx/km

print(tc," sec  ", rmx*100, " cm")


    
#plt.plot(100*R,V)
#plt.xlabel('Radius in cm')
#plt.ylabel('Velocity in m/s')
#plt.grid()
#plt.show    


time = np.linspace(0,1*t0(850,rmx),500)
velocity_vt = []
pos = []
for t in time:
    velocity = f_v(t,rmx)
    x = f_x(t,rmx)
    
    velocity_vt.append(velocity)
    pos.append(x)

plt.plot(time,velocity_vt)
plt.xlabel('Time in s')
#plt.xlabel('Distance in m')
plt.ylabel('Velocity in m/s')
plt.grid()
plt.show
   









#Rnew = rmx/1.5
#time = np.linspace(0,1*t0(850,Rnew),500)
#velocity_vt = []
#pos = []
#for t in time:
#    velocity = f_v(t,Rnew)
#    x = f_x(t,Rnew)
#    
#    velocity_vt.append(velocity)
#    pos.append(x)
#
#plt.plot(time,velocity_vt)
#plt.xlabel('Time in s')
##plt.xlabel('Distance in m')
#plt.ylabel('Velocity in m/s')
#plt.grid()
#plt.show








