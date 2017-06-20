import numpy as np
import dynamics2 as d2
#import forceTorque as ft

def rk4(f,t0,x0,h):
#range kutta order 4 solver
	k1 = h*f(t0,x0)
	k2 = h*(ft0+h/2,x0+k1/2)
	k3 = h*f(t0+h/2,x0+k2/2)
	k4 = h*f(t0+h,x0+k3)
	x1 = x0.copy
	x1 += (k1 + 2*k2 + 2*k3 + k4)/6

	return x1

def rk42(x0,h,F,T):
#range kutta order 4 solver, assuming magnetic field does not change in time t to t+h
	k1 = h*d2.dynamics2(x0,F,T)
	k2 = h*d2.dynamics2(x0+k1/2,F,T)
	k3 = h*d2.dynamics2(x0+k2/2,F,T)
	k4 = h*d2.dynamics2(x0+k3,F,T)
	x1 = x0.copy
	x1 += (k1 + 2*k2 + 2*k3 + k4)/6

	return x1

def simpsonG(F,T,k,g_dFdT,dm,i,v_dL,v_pos_sat_b,q):
	
	v_pos_dL_b = i*v_dL
	v_dF_i, v_dT_b = k*g_dFdT(dm,v_pos_sat_b,v_pos_dL_b,q)
	F1 = F + v_dF_i
	T1 = T + v_dT_b

	return F1,T1

def simpsonM(F,T,k,m_dFdT,i,v_pos_sat_b,v_dL,v_dL_cap_i,v_v_sat_i,q,t,dL):

	v_pos_dL_b = v_pos_sat_b + i*v_dL
	v_dF_i, v_dT_i = k*m_dFdT(v_pos_dL_b,v_dL_cap_i,v_v_sat_i,q,t,dL)
	F1 = F + v_dF_i
	T1 = T + v_dT_b

	return F1, T1


