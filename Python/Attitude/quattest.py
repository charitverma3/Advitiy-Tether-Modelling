import qnv
import numpy as np
from math import sin,cos,sqrt,pi
from constants import G,M,R
def dqdt(q,w):
	q_dot = 0.5*qnv.quatMultiply(q,w)
	q_dot = q_dot.reshape((4,1))

	return q_dot

def drdt(p):
	v_dot = -G*M*p[0:3].reshape(3,1)/np.linalg.norm(p[0:3])**3
	return np.vstack((p[3:6].reshape(3,1),v_dot))

def drdt2(p,a):
	v_dot = a
	return np.vstack((p[3:6].reshape(3,1),v_dot))

def rk4(x0,h,w):
#range kutta order 4 solver
	k1 = h*dqdt(x0,w)
	k2 = h*dqdt(x0+k1/2,w)
	k3 = h*dqdt(x0+k2/2,w)
	k4 = h*dqdt(x0+k3,w)
	x1 = x0.copy()
	#print x1, k1
	x1 = x1 + (k1 + 2*k2 + 2*k3 + k4)/6

	x1 = x1/np.linalg.norm(x1)
	return x1
def rk4a(p,h):
	k1 = h*drdt(p)
	k2 = h*drdt(p+k1/2)
	k3 = h*drdt(p+k2/2)
	k4 = h*drdt(p+k3)
	x1 = p.copy()
	#print x1, k1
	x1 = x1 + (k1 + 2*k2 + 2*k3 + k4)/6

	return x1

def rk4b(p,h):
	a = -G*M*p[0:3].reshape(3,1)/np.linalg.norm(p[0:3])**3
	k1 = h*drdt2(p,a)
	k2 = h*drdt2(p+k1/2,a)
	k3 = h*drdt2(p+k2/2,a)
	k4 = h*drdt2(p+k3,a)
	x1 = p.copy()
	#print x1, k1
	x1 = x1 + (k1 + 2*k2 + 2*k3 + k4)/6

	return x1

theta0 = 0
pos0 = np.array([[R + 600e3],[0.],[0.]])
dist0 = np.linalg.norm(pos0)
v0 = np.array([[0],[sqrt(G*M/dist0)],[0]])
w0 = np.array([[0],[-1*sqrt(G*M/dist0**3)],[0]])
w0n = np.linalg.norm(w0)
v_L_b = np.array([[0.], [0.], [1.]])
m_eu0 = np.array([[0,0,-1], [1,0,0], [0,-1,0]])
q0 = qnv.rotm2quat(m_eu0)
q0 = q0.reshape((4,1))

time_i = 0
time_f = 10*pi/w0n
step_size = 0.1
nT = int((time_f - time_i)/step_size)

dot = np.zeros(nT+1)
dot[0] = np.dot((qnv.quatRotate(q0,v_L_b)).T, pos0)/(np.linalg.norm(pos0))
r = np.zeros(nT+1)
r[0] = np.linalg.norm(pos0)
q1 = q0
theta1 = theta0
pos1 = pos0
s1 = np.vstack((pos0,v0))

for n in range(0,nT):
	if n%10000 == 0:
		print n*step_size
	q2 = rk4(q1,step_size,np.vstack(([0.],w0)))
	s2 = rk4a(s1, step_size)
	theta2 = theta1 + w0n*step_size
	#pos2 = np.array([[cos(theta2)],[sin(theta2)],[0.]]) 
	pos2 = s2[0:3].reshape(3,1)
	r[n+1] = np.linalg.norm(pos2)
	#print pos2
	dot[n+1] = np.dot((qnv.quatRotate(q2,v_L_b)).T, pos2)/(np.linalg.norm(pos2))
	q1 = q2
	theta1 = theta2
	s1 = s2