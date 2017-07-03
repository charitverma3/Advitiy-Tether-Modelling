import numpy as np
import qnv as qnv
import frames as fs
from constants import *
from pyigrf12 import runigrf12
from solver import *

def gravityForceTorque(state):
    #gives force in ECIF and torque in body frame due to gravity
    v_dL = v_L_b/nLg
    dL = np.linalg.norm(v_dL)
    dm = mu_m*dL
    v_F_i = np.zeros((3,1))
    v_T_b = np.zeros((3,1))
    q = state[6:10].copy()
    qi = qnv.quatInv(q);
    v_pos_sat_i = state[0:3].copy()
    v_pos_sat_b = qnv.quatRotate(qi,v_pos_sat_i)
    
    if nLg>1 :
        for i in range(1,int(nLg/2+1)) :

            v_F_i, v_T_b = simpsonG(v_F_i,v_T_b,4,g_dFdT,dm,(2*i-1),v_dL,v_pos_sat_b,q)
            
        for i in range(2,int(nLg/2+1)) : 

            v_F_i, v_T_b = simpsonG(v_F_i, v_T_b,2,g_dFdT,dm,(2*i-2),v_dL,v_pos_sat_b,q)
            
        for i in [1,int(nLg)]:
            
            v_F_i, v_T_b = simpsonG(v_F_i, v_T_b,1,g_dFdT,dm,i,v_dL,v_pos_sat_b,q) 

        v_F_i = v_F_i/3
        v_T_b = v_T_b/3

    else:
        v_F_i, v_T_b = simpsonG(v_F_i, v_T_b,1,g_dFdT,dm,nLg/2,v_dL,v_pos_sat_b,q)  
        

    return v_F_i, v_T_b

def g_da(p1,p2):
    #gives gravitational acceleration at position p1+p2
    da = -G*M*(p1+p2)/(np.linalg.norm(p1+p2)**3)
    return da

def g_dFdT(dm,v_pos_sat_b,v_pos_dL_b,q):
    #gives differential force in ECIF and differential torque in body frame due to gravity
    v_dF_b = dm*(g_da(v_pos_sat_b,v_pos_dL_b))
    v_dF_i = qnv.quatRotate(q,v_dF_b)
    v_dT_b = qnv.cross1(v_pos_dL_b,v_dF_b)
    return v_dF_i, v_dT_b


def magneticForceTorque(state,t):
    q = state[6:10].copy()
    qi = qnv.quatInv(q)
    v_pos_sat_i = state[0:3].copy()
    v_pos_sat_b = qnv.quatRotate(qi, v_pos_sat_i)
    v_v_sat_i = state[3:6].copy()
    v_L_i = qnv.quatRotate(q,v_L_b)
    v_dL_i = v_L_i/nLb
    v_dL_cap_i = v_dL_i/np.linalg.norm(v_dL_i)
    dL = np.linalg.norm(v_dL_i)
    v_dL_b = v_L_b/nLb
    v_F_i = np.array([[0.], [0.], [0.]])
    v_T_b = np.array([[0.], [0.], [0.]])

    if nLb>1 :
        for i in range(1,nLb/2+1):
            v_F_i, v_T_b = simpsonM(v_F_i, v_T_b,4,m_dFdT,(2*i-1),v_pos_sat_b,v_dL_b,v_dL_cap_i,v_v_sat_i,q,t,dL) 
            
        for i in range(2,nLb/2+1):
            v_F_i, v_T_b = simpsonM(v_F_i, v_T_b,2,m_dFdT,(2*i-2),v_pos_sat_b,v_dL_b,v_dL_cap_i,v_v_sat_i,q,t,dL)
            
        for i in [1,nLg]:
            v_F_i, v_T_b = simpsonM(v_F_i, v_T_b,1,m_dFdT,i,v_pos_sat_b,v_dL_b,v_dL_cap_i,v_v_sat_i,q,t,dL)

        v_F_i = v_F_i/3
        v_T_b = v_T_b/3

    else:
        v_F_i, v_T_b = simpsonM(v_F_i, v_T_b,1,m_dFdT,nLb/2,v_pos_sat_b,v_dL_b,v_dL_cap_i,v_v_sat_i,q,t,dL)

    return v_F_i, v_T_b


def m_dFdT(v_pos_dL_b,v_dL_cap_i,v_v_sat_i,q,t,dL):
    qi = qnv.quatInv(q)
    v_pos_dL_i = qnv.quatRotate(q,v_pos_dL_b)
    height = np.linalg.norm(v_pos_dL_i) - R
    height = height/1e3
    v_pos_dL_e = fs.ecif2ecef(v_pos_dL_i,t)
    lat, lon = fs.latlon(v_pos_dL_e.reshape((3,)))
    B = runigrf12(day,0,1,height,lat,lon)
    v_B_n = np.vstack((B[0],B[1], B[2]))
    v_B_n = v_B_n*1e-9 #convert from nT to T
    v_B_e = fs.ned2ecef(v_B_n.reshape((3,)),lat,lon)
    v_dL_cap_e = fs.ecif2ecef(v_dL_cap_i,t)
    v_v_sat_e = fs.ecif2ecef(v_v_sat_i,t)
    e = qnv.dot1(qnv.cross1(v_v_sat_e, v_B_e), v_dL_cap_e)
    i = e/mu_r
    v_dF_e = dL*i*qnv.cross1(v_dL_cap_e,v_B_e)
    v_dF_i = fs.ecef2ecif(v_dF_e, t)
    v_dF_b = qnv.quatRotate(qi, v_dF_i)
    v_dL_cap_b = v_L_b/np.linalg.norm(v_L_b)
    v_dT_b = qnv.cross1(v_dL_cap_b, v_dF_b)
    return v_dF_i, dL*v_dT_b


def getB(state,t):
    q = state[6:10].copy()
    qi = qnv.quatInv(q)
    v_pos_sat_i = state[0:3].copy()
    v_pos_sat_b = qnv.quatRotate(qi, v_pos_sat_i)
    v_pos_dL_i = qnv.quatRotate(q,v_pos_sat_b + v_L_b/2.0)
    height = np.linalg.norm(v_pos_dL_i) - R
    height = height/1e3
    v_pos_dL_e = fs.ecif2ecef(v_pos_dL_i,t)
    lat, lon = fs.latlon(v_pos_dL_e.reshape((3,)))
    B = runigrf12(day,0,1,height,lat,lon)
    v_B_n = np.vstack((B[0],B[1], B[2]))
    v_B_n = v_B_n*1e-9 #convert from nT to T
    v_B_e = fs.ned2ecef(v_B_n.reshape((3,)),lat,lon)
    return v_B_e


def m_dFdT2(v_pos_dL_b,v_dL_cap_i,v_v_sat_i,q,t,dL,v_B_e):
    qi = qnv.quatInv(q)
    v_pos_dL_i = qnv.quatRotate(q,v_pos_dL_b)
    v_pos_dL_e = fs.ecif2ecef(v_pos_dL_i,t)
    v_dL_cap_e = fs.ecif2ecef(v_dL_cap_i,t)
    v_v_sat_e = fs.ecif2ecef(v_v_sat_i,t)
    e = qnv.dot1(qnv.cross1(v_v_sat_e, v_B_e), v_dL_cap_e)
    i = e/mu_r
    v_dF_e = dL*i*qnv.cross1(v_dL_cap_e,v_B_e)
    v_dF_i = fs.ecef2ecif(v_dF_e, t)
    v_dF_b = qnv.quatRotate(qi, v_dF_i)
    v_dL_cap_b = v_L_b/np.linalg.norm(v_L_b)
    v_dT_b = qnv.cross1(v_dL_cap_b, v_dF_b)
    return v_dF_i, dL*v_dT_b

def magneticForceTorque2(state,t,B):
    q = state[6:10].copy()
    qi = qnv.quatInv(q)
    v_pos_sat_i = state[0:3].copy()
    v_pos_sat_b = qnv.quatRotate(qi, v_pos_sat_i)
    v_v_sat_i = state[3:6].copy()
    v_L_i = qnv.quatRotate(q,v_L_b)
    v_dL_i = v_L_i/nLb
    v_dL_cap_i = v_dL_i/np.linalg.norm(v_dL_i)
    dL = np.linalg.norm(v_dL_i)
    v_dL_b = v_L_b/nLb
    v_F_i = np.array([[0.], [0.], [0.]])
    v_T_b = np.array([[0.], [0.], [0.]])

    if nLb>1 :
        for i in range(1,nLb/2+1):
            v_F_i, v_T_b = simpsonM2(v_F_i, v_T_b,4,m_dFdT2,(2*i-1),v_pos_sat_b,v_dL_b,v_dL_cap_i,v_v_sat_i,q,t,dL,B) 
            
        for i in range(2,nLb/2+1):
            v_F_i, v_T_b = simpsonM2(v_F_i, v_T_b,2,m_dFdT2,(2*i-2),v_pos_sat_b,v_dL_b,v_dL_cap_i,v_v_sat_i,q,t,dL,B)
            
        for i in [1,nLg]:
            v_F_i, v_T_b = simpsonM2(v_F_i, v_T_b,1,m_dFdT2,i,v_pos_sat_b,v_dL_b,v_dL_cap_i,v_v_sat_i,q,t,dL,B)

        v_F_i = v_F_i/3
        v_T_b = v_T_b/3

    else:
        v_F_i, v_T_b = simpsonM2(v_F_i, v_T_b,1,m_dFdT2,nLb/2,v_pos_sat_b,v_dL_b,v_dL_cap_i,v_v_sat_i,q,t,dL,B)

    return v_F_i, v_T_b













