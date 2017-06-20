import numpy as np
import forceTorque as ft
import solver as slv
#simulation variables

time_i = 0
time_f = 1
step_size = 0.1
nT = (time_f - time_i)/step_size
state = np.zeros((13,nT+1))
state(:,0) = np.array([])
time = time_i

for n in (1:nT+1):
	#Fg, Tg = ft.gravityForceTorque(state)
	Fm, Tm = ft.magneticForceTorque(state,time)
	state(:,n+1) = slv.rk42(state(:,n), step_size, Fm, Tm)






