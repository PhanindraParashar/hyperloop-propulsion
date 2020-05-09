import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# function that returns dy/dt
def model(theta,t):
    
    m = .2
    mu = .5
    g = 10
    l = 15/100
    r = 7/100
    
    gratio = 2.5
    
    Tst = 3.5/gratio
    Wrpm = 6500*gratio
    
    Wmx = Wrpm*(np.pi/30)
    km = (Tst/Wmx)
    
    beta = np.arcsin(np.sin(theta)*(r/l))
    alpha = np.pi - (beta + theta)
    
    Fp = 34
    
    A = Fp*r*np.sin(theta) - mu*m*g*r*np.sin(theta)
    B = np.cos(beta) + mu*np.sin(beta)
    C = m*r*np.sin(theta)
    D = r*np.cos(theta) + (r*r/l)*np.cos(2*theta)
    E = C*D
    F = A - B*Tst
    
    
    det = np.sqrt((B*km)**2 - 4*E*F)
    
    dydt = (det - (B*km))/(2*E)
    
    return dydt

# initial condition
theta0 = (np.pi/100)

# time points
t = np.linspace(0,.04,400)

# solve ODE
theta = odeint(model,theta0,t)

# plot results
plt.plot(t,theta)
plt.xlabel('time')
plt.ylabel('theta(t)')
plt.show()