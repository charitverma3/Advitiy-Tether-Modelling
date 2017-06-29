import numpy as np

def dynamics(state,t):
	#dynamic equations of system
	v_Fg_i, v_Tg_b = gravityForceTorque(state)
	v_Fb_i, v_Tb_b = magneticForceTorque(state)

	v_F_i = v_Fb_i + v_Fg_i
	v_T_b = v_Tb_b + v_Tg_b

	a = v_F_i/Ms #acceleration in inertial frame
	print Ms
	q = state[7:11].copy()
	omega = state[11:14].copy()

	q_dot = 0.5*q*omega
	omega_dot = inv(m_Inertia)*(v_T_b - cross1(omega, m_Inertia*omega))

	y = np.vstack(([state[4:7],a,q_dot,omega_dot])

	return y