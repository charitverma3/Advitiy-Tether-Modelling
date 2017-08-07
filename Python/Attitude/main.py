import numpy as np
import forceTorque as ft
import solver as slv
import math
import qnv
from constants import * 
import scipy.io as sio
import time
import datetime
import os
import satellite
#simulation variables

dir_now = os.path.normpath(datetime.datetime.now().strftime('%Y-%m-%d %H-%M-%S'))

t1 = time.time()
time_i = 0.
#time_f = math.pi/(2*math.sqrt(G*M/R**3))
time_f = 10.
step_size = 0.1
nT = int((time_f - time_i)/step_size)
s_time = time_i

state = np.zeros((13,nT+1)) #state = (pos from earth in ECIF, velocity, quaternion, angular velocity wrt ECIF in body frame) quaternion rotates body frame vector into inertial frame and defined as (scalar,vector)
dot = np.zeros(nT+1)
energy = np.zeros(nT+1)
r = np.zeros((nT+1))
omega = np.zeros((nT+1))

state[:,0] = state0.reshape((1,13))
r[0] = np.linalg.norm(state0[0:3])
dot[0] = np.dot((qnv.quatRotate(q0,v_L_b)).T, state0[0:3])/(L*np.linalg.norm(state0[0:3]))
energy[0] = 0.5*Ms*(np.linalg.norm(state0[3:6]))**2 - G*M*Ms/(np.linalg.norm(state0[0:3])) 
omega[0] = np.linalg.norm(w0)

sat = satellite.Satellite(state0,time_i)
print sat

for n in range(0,nT):
	if (n%10000==0):
		t2 = time.time()
		t3 = t2 - t1
		print n*step_size
		print t3
		print dot[n]
	#Fg, Tg = ft.gravityForceTorque(state)
	state_now = sat.getState()
	B = ft.getB(sat)
	print "x"
	#Fm, Tm = ft.magneticForceTorque(state_now,s_time)
	#Fm = np.zeros((3,1))
	#Tm = np.zeros((3,1))

	state_next = (slv.rk43(sat, step_size, B)).reshape((13,1))
	s_time = s_time + step_size
	sat.setState(state_next)
	sat.setTime(s_time)
	r[n+1] = sat.getR()
	dot[n+1] = np.dot((sat.v_dL_cap_i, sat.getPos()))/(r[n+1])
	energy[n+1] = sat.getEnergy()
	omega[n+1] = np.linalg.norm(sat.getW())

t2 = time.time()
t3 = t2 - t1
print t3


os.chdir('Logs-With-Torque-2')
os.mkdir(dir_now)
os.chdir(dir_now)
sio.savemat('state.mat', mdict={'state':state})
sio.savemat('r.mat', mdict={'r':r})
sio.savemat('dot.mat', mdict={'dot':dot})
sio.savemat('energy.mat', mdict={'energy':energy})
sio.savemat('omega.mat', mdict={'omega':omega})


