from constants import *
import numpy as np
import qnv
import frames as fs

class Satellite:

	def __init__(self,state0,time0):

		self.time = time0
		self.setState(state0)		
		print "init"
	def setState(self,state):

		self.state = state.copy()
		v_pos_sat_i = state[0:3].copy()
		v_v_sat_i = state[3:6].copy()
		q = state[6:10],copy()
		self.qi = qnv.quatInv(q)
		self.v_pos_sat_b = qnv.quatRotate(qi, v_pos_sat_i)
		v_L_i = qnv.quatRotate(q,v_L_b)
		self.v_dL_cap_i = v_L_i/np.linalg.norm(v_L_i)
		self.v_dL_cap_e = fs.ecif2ecef(v_dL_cap_i,self.time)
		self.v_v_sat_e = fs.ecif2ecef(v_v_sat_i,self.time)

	def getState(self):

		return self.state

	def setPos(self,pos):

		self.state[0:3] = pos

	def getPos(self):

		return self.state[0:3]

	def setVel(self,v):

		self.state[3:6] = v

	def getVel(self):

		return self.state[3:6]

	def setQ(self,q):

		self.state[6:10] = q

	def getQ(self):

		return self.state[6:10]

	def getQi(self):

		return self.qi

	def setW(self,omega):

		self.state[10:13] = omega

	def getW(self):

		return self.state[10:13]

	def getEnergy(self):
		pos = self.getPos()
		v = self.getVel()
		omega = self.getW()
		T = 0.5*Ms*(np.linalg.norm(v))**2 - G*M*(Ms + mu_m*L/(np.linalg.norm(pos)) + 0.5*np.dot(omega, np.matmul(m_Inertia,omega))

	def setEmf(self,e):

		self.emf = e

	def getEmf(self):
		
		return self.emf

	def setTime(self,y):
		self.time = t

	def getTime(self):
		return self.time