import numpy as np
import datetime as dt
import math
import qnv
w_earth = 7.2921159e-5; #rad per second
G=6.67e-11; #universal gravitational constant, SI
M=5.972e24; #mass of earth, kg
R=6.3718e6; #radius of earth, m
day = dt.datetime(2017,5,30)
L = 100. #length of tether in m
v_L_b = L*np.array([[0.], [0.], [1.]])
b = 1e-3 #width of tether
mu_m = 0*0.27/L #mass per unit length of tether
nLb = 1.
nLg = 1.
mu_r = 0.1
Ms = 10.

Ixxs = .17007470856
Iyys = .17159934710
Izzs = .15858572070
Ixys = .00071033134
Iyzs = .00240388659
Ixzs = .00059844292

Ixxt = mu_m*L*(b**2/12 + L**2/3)
Iyyt = Ixxt
Izzt = mu_m*L*(b**2/6)
Ixyt = 0.
Ixzt = 0.
Iyzt = 0.

Ixx = Ixxs + Ixxt
Iyy = Iyys + Iyyt
Izz = Izzs + Izzt
Ixy = Ixys + Ixyt
Ixz = Ixzs + Ixzt
Iyz = Iyzs + Iyzt

m_Inertia = np.array([[Ixx, -1*Ixy, -1*Ixz], [-1*Ixy, Iyy, -1*Iyz], [-1*Ixz, -1*Iyz, Izz]])
m_Inertia_inv = np.linalg.inv(m_Inertia)

incl = math.radians(98)

m_eu98 = np.array([[0.,0.,-1.], [math.cos(incl),math.sin(incl),0.], [math.sin(incl),-math.cos(incl),0.]])
#m_eu0 = np.array([[0.,0.,-1.],[1.,0.,0.],[0,-1.,0.]])
q0 = qnv.rotm2quat(m_eu98)
q0 = q0.reshape((4,1))
pos0 = np.array([[R + 6e3],[0.],[0.]])
dist0 = np.linalg.norm(pos0)
#v0 = np.array([[0.],[math.sqrt(G*M/dist0)],[0.]])
v0 = np.array([[0],[math.sqrt(G*M/dist0)*math.cos(incl)],[math.sqrt(G*M/dist0)*math.sin(incl)]])
w0 = np.array([[0.],[-1*math.sqrt(G*M/dist0**3)],[0.]]) 
#w0 = np.array([[0.], [0.], [0.]])
state0 = np.vstack((pos0,v0,q0,w0))







