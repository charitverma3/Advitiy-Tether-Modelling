import numpy as np
import qnv
from constants import *
import forceTorque as fT
def dynamics3(sat,B):
#dynamic equations when magnetic field is assumend constant between t and t+dt for computational ease 
	v_Fg_i, v_Tg_b = fT.gravityForceTorque(sat)
	v_Fs_i = Ms*fT.g_da(sat.getPos(),0)
	v_Fb_i, v_Tb_b = fT.magneticForceTorque(sat,B)
	v_F_i = v_Fb_i + v_Fg_i + v_Fs_i
	v_T_b = v_Tb_b + v_Tg_b

	a = v_F_i/(Ms + mu_m*L) #acceleration in inertial frame
	#print np.dot(a,state[0:3])
	q = sat.getQ()
	omega = sat.getW()
	omega_q = np.hstack(([0.],omega))
	q_dot = 0.5*qnv.quatMultiply(q,omega_q)
	omega_dot = np.dot(m_Inertia_inv,(v_T_b - qnv.cross1(omega, np.dot(m_Inertia,omega))))
	omega_dot = omega_dot*0.
	y = np.hstack((sat.getVel(),a,q_dot,omega_dot))
	
	return y