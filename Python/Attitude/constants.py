import numpy as np
import datetime as dt
import math
import qnv
w_earth = 7.2921159e-5; #rad per second
G=6.67e-11; #universal gravitational constant, SI
M=5.972e24; #mass of earth, kg
R=6371.8e3; #radius of earth, m
day = dt.datetime(2017,5,30)
L = 100. #length of tether in m
v_L_b = L*np.array([[0.], [0.], [1.]])

Ixx = .17007470856
Iyy = .17159934710
Izz = .15858572070
Ixy = .00071033134
Iyz = .00240388659
Ixz = .00059844292
m_Inertia = np.array([[Ixx, -1*Ixy, -1*Ixz], [-1*Ixy, Iyy, -1*Iyz], [-1*Ixz, -1*Iyz, Izz]])
m_Inertia_inv = np.linalg.inv(m_Inertia)


nLb = 1.
nLg = 1.
Ms = 10.
mu_m = 0.27/L
mu_r = 0.1

incl = math.radians(98)

m_eu98 = np.array([[0.,0.,-1.], [math.cos(incl),math.sin(incl),0.], [math.sin(incl),-math.cos(incl),0.]])
#m_eu0 = np.array([[0.,0.,-1.],[1.,0.,0.],[0,-1.,0.]])
q0 = qnv.rotm2quat(m_eu98)
q0 = q0.reshape((4,1))
pos0 = np.array([[R + 5.1e5],[0.],[0.]])
dist0 = np.linalg.norm(pos0)
#v0 = np.array([[0.],[math.sqrt(G*M/dist0)],[0.]])
v0 = np.array([[0],[math.sqrt(G*M/dist0)*math.cos(math.radians(98))],[math.sqrt(G*M/dist0)*math.sin(math.radians(98))]])
w0 = np.array([[0.],[-1*math.sqrt(G*M/dist0**3)],[0.]]) 
#w0 = np.array([[0.], [0.], [0.]])
state0 = np.vstack((pos0,v0,q0,w0))







