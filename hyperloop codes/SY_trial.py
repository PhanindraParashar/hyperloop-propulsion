import numpy as np

theta = np.pi/4

m = .05
mu = .4
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
    
A = np.sin(alpha) + mu*np.cos(alpha)
B = m*r*mu*g
D = A*Tst - B
C = (m*(r**3)/l)*np.cos(2*theta) + (m*r*r*np.cos(theta))
    
det = np.sqrt((A*km)**2 - 4*C*D)
    
dydt = (det + (A*km))/(2*C)

print(dydt)