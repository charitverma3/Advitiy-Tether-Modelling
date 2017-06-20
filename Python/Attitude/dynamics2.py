import numpy as np
import qnv
from constants import *

def dynamics(state,v_Fb_i, v_Tb_b):
#dynamic equations when magnetic field is assumend constant between t and t+dt for computational ease 
	v_Fg_i, v_Tg_b = gravityForceTorque(state)
	
	v_F_i = v_Fb_i + v_Fg_i
	v_T_b = v_Tb_b + v_Tg_b

	a = v_F_i/Ms #acceleration in inertial frame

	q = state[7:11]
	omega = state[11:14]
	omega_q = np.vstack(([0.],omega))
	q_dot = 0.5*qnv.quatMultiply(q,omega_q)
	omega_dot = np.dot(m_Inertia_inv,(v_T_b - qnv.cross1(omega, np.dot(m_Inertia,omega))))

	y = np.vstack((state[4:7],a,q_dot,omega_dot))

	return y