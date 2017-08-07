import numpy as np
import dynamics2 as d2
import dynamics3 as d3
#import forceTorque as ft
import copy
def rk4(f,t0,x0,h):
#range kutta order 4 solver
	k1 = h*f(t0,x0)
	k2 = h*(ft0+h/2,x0+k1/2)
	k3 = h*f(t0+h/2,x0+k2/2)
	k4 = h*f(t0+h,x0+k3)
	x1 = x0.copy()
	x1 = x1 + (k1 + 2*k2 + 2*k3 + k4)/6

	return x1

def rk42(x0,h,F,T):
#range kutta order 4 solver, assuming magnetic force does not change in time t to t+h
	k1 = h*d2.dynamics2(x0,F,T)
	k2 = h*d2.dynamics2(x0+k1/2,F,T)
	k3 = h*d2.dynamics2(x0+k2/2,F,T)
	k4 = h*d2.dynamics2(x0+k3,F,T)
	x1 = x0.copy()
	#print k3
	x1 = x1 + (k1 + 2*k2 + 2*k3 + k4)/6

	x1[6:10] = x1[6:10]/np.linalg.norm(x1[6:10])
	return x1

def rk43(sat,h,B):
#range kutta order 4 solver, assuming magnetic field in ECEF does not change in time t to t+h
	s = copy.deepcopy(sat)
	x0 = s.getState()
	k1 = h*d3.dynamics3(s,B)
	s.setState(x0+k1/2)
	k2 = h*d3.dynamics3(s,B)
	s.setState(x0+k2/2)
	k3 = h*d3.dynamics3(s,B)
	s.setState(x0+k3)
	k4 = h*d3.dynamics3(s,B)

	x1 = x0 + (k1 + 2*k2 + 2*k3 + k4)/6

	x1[6:10] = x1[6:10]/np.linalg.norm(x1[6:10])
	return x1

def simpsonG(F,T,k,g_dFdT,i,sat):
	
	
	v_dF_i, v_dT_b = g_dFdT(sat,i)
	F1 = F + k*v_dF_i
	T1 = T + k*v_dT_b

	return F1,T1

def simpsonM(F,T,k,m_dFdT,i,sat,B):

	
	v_dF_i, v_dT_b = m_dFdT(sat,B)
	F1 = F + k*v_dF_i
	T1 = T + k*v_dT_b

	return F1, T1

def simpsonM2(F,T,k,m_dFdT2,i,v_pos_sat_b,v_dL,v_dL_cap_i,v_v_sat_i,q,t,dL,B):

	v_pos_dL_b = v_pos_sat_b + i*v_dL
	v_dF_i, v_dT_b = m_dFdT2(v_pos_dL_b,v_dL_cap_i,v_v_sat_i,q,t,dL,B)
	F1 = F + k*v_dF_i
	T1 = T + k*v_dT_b

	return F1, T1


