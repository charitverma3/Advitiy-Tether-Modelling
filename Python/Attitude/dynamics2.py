import numpy as np

def dynamics(state,v_Fb_i, v_Tb_b):
#dynamic equations when magnetic field is assumend constant between t and t+dt for computational ease 
	v_Fg_i, v_Tg_b = gravityForceTorque(state)
	
	v_F_i = v_Fb_i + v_Fg_i
	v_T_b = v_Tb_b _ v_Tg_b

	a = v_F_i/Ms #acceleration in inertial frame

	q = state[7:11]
	omega = state[11:13]

	q_dot = 0.5*q*omega
	omega_dot = inv(m_Inertia)*(v_T_b - cross(omega, m_Inertia*omega))

	y = [state(4:6);a;q_dot;omega_dot]

	return y