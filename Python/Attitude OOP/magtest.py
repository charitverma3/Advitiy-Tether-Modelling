from constants import G,M,R,L
import forceTorque as ft
import frames as fs
import qnv
import numpy as np
from math import sqrt
pos0 = np.array([[R],[0],[0]])
dist0 = np.linalg.norm(pos0)
#v0 = np.array([[0],[sqrt(G*M/dist0)],[0]])
v0 = np.array([[3e3],[3e3],[0]])

w0 = np.array([[0],[-1*sqrt(G*M/dist0**3)],[0]])
w0n = np.linalg.norm(w0)

m_eu0 = np.array([[0,0,-1], [1,0,0], [0,-1,0]])
q0 = qnv.rotm2quat(m_eu0)
q0 = q0.reshape((4,1))
qi = qnv.quatInv(q0)

v_L_b = L*np.array([[0.], [0.], [1.]])
v_L_i = qnv.quatRotate(q0,v_L_b)

t=100

pos1 = qnv.quatRotate(qi,pos0)
v_pos_dL_b = pos1 + v_L_b
v_i = v0
dL = L

a,b = ft.m_dFdT(v_pos_dL_b,v_L_i/100,v_i,q0,t,dL)
a1 = qnv.quatRotate(qi,a)
print a1,b