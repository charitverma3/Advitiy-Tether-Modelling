import numpy as np
import datetime as dt

w_earth = 7.2921159e-5; #rad per second
G=6.67e-11; #universal gravitational constant, SI
M=5.972e24; #mass of earth, kg
R=6371.8e3; #radius of earth, m
day = dt.datetime(2017,6,17)
L = 100 #length of tether in m
v_L_b = L*np.array([[0.], [0.], [1.]])

Ixx = .17007470856
Iyy = .17159934710
Izz = .15858572070
Ixy = .00071033134
Iyz = .00240388659
Ixz = .00059844292
m_Inertia = np.array([[Ixx, -1*Ixy, -1*Ixz], [-1*Ixy, Iyy, -1*Iyz], [-1*Ixz, -1*Iyz, Izz]])
m_Inertia_inv = np.linalg.inv(m_Inertia)


nLb = 1
nLg = 1
pi = np.pi
Ms = 10
mu_m = 0.01
mu_r = 0.1









