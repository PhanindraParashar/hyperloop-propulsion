import numpy as np
import matplotlib.pyplot as plt


def f_x(t,r,m) :
    Fmax = Tst/r
    f = mu*m*g
    Vlmax = r*Wmax
    Ce = Fmax/(m*Vlmax)
    
    eta = (1-(f/Fmax))
    return eta*Vlmax*t - (eta*m/Fmax)*(Vlmax**2)*(1-np.exp(-Ce*t))

def f_v(t,r,m) :
    Fmax = Tst/r
    f = mu*m*g
    Vlmax = r*Wmax
    Ce = Fmax/(m*Vlmax)
    
    eta = (1-(f/Fmax))
    return eta*Vlmax*(1-np.exp(-Ce*t))


def t0(s,r,m):
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

gearR = .11
wmax_rpm = 17900*gearR
Tst = 639/gearR  # Nm

mu = 0.015
g = 9.8


Wmax = wmax_rpm*(np.pi/30)

km = Tst/Wmax
rcm = 5.75
r = rcm/100

M = np.linspace(50,60,1000)  # kg
V = []
for m in M:
    t = t0(850,r,m)
    vel = f_v(t,r,m)
    V.append(vel)

#plt.plot(M,V)
#plt.xlabel("Mass in Kg")
#plt.ylabel("Velocity in m/s")
#plt.grid(True)



#af = 0.0035*1
af = .34

m = 1847
r = .512

A = (Tst/(m*r))-mu*g

km = Tst/Wmax

B = km/(m*r*r)

C = af/m

D = np.sqrt(B**2 + 4*A*C)

def v_ad(t):
    return (D*np.tanh(D*t/2 + np.arctanh(B/D)) - B)/(2*C) 

def x_ad(t):
    Th = np.tanh(D*t/2)
    
    t1 = (1/C)*np.log(1+(B*Th/D))
    t2 = (-2*A/D)*(1/(B+D))*np.log(1-Th)
    t3 = (-2*A/D)*(1/(D-B))*np.log(1+Th)
    return t1+t2+t3

time = np.linspace(3,5,1000)
ve = []
for t in time:
    vel = v_ad(t)
    ve.append(vel)
plt.plot(time,ve)
plt.grid(True)
