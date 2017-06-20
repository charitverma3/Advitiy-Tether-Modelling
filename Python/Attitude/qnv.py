import numpy as np

def cross1(v1,v2):
	#for cross product of 2 column vectors

	v3 = v1.reshape((1,3))
	v4 = v2.reshape((1,3))
	v5 = np.cross(v3,v4)
	v6 = v5.reshape((3,1))
	return v6


def quatInv(q):
	#to get inverse of a quaternion
	qi = np.vstack((q[0:1],-1*q[1:4]))
	qi = qi/np.linalg.norm(q)
	return qi


def quatMultiply(q1,q2):

	#quaternion is scalar, vector. function multiplies 2 quaternions

	a1 = q1[0:1].copy()
	a2 = q2[0:1].copy()

	b1 = (q1[1:4].copy()).reshape(3,)
	b2 = (q2[1:4].copy()).reshape(3,)
	
	a = a1 + a2 - np.dot(b1,b2)
	b = a1*b2 + a2*b1 + np.cross(b1,b2)

	b = b.reshape((3,1))

	q = np.vstack((a,b))

	return q

def quatRotate(q,x):
	
	#rotates vecctor x by quaternion q
	qi = quatInv(q)
	y = np.vstack(([0.],x.copy()))
	y = quatMultiply(q,y)
	y = quatMultiply(y,qi)

	x2 = y[1:4]
	return x2

if __name__ == "__main__":
	pass



