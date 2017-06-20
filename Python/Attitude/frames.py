import numpy as np
from constants import *

def latlon(x):
	#get the latitude and longitude given position in ECEF
	lat=sgn(x[2])*(acos(((x[0]^2+x[1]^2)^0.5)/((x[0]^2+x[1]^2+x[2]^2)^0.5)))*90/(np.pi/2) 
        
	# longitude calculation given position, lon is longitude
	if x[1]==0:
		if x[0]>=0 :
			lon = 0
		else:
			lon = np.pi

	else:
		lon=sgn(x[1])*acos(x[0]/((x[0]^2+x[1]^2)^0.5))*90/(np.pi/2)     #x,y,z>0 

	#x axis is intersection of 0 longitude and 0 latitutde

def sgn(x):
	#signum function
	if(x==0):
		y = 0
	elif (x>0):
		y = 1
	else
		y = -1

	return y

def ecif2ecef(x,t):
	theta = w_earth*t #in radian
	DCM = np.array([[cos(theta), sin(theta), 0], [-1*sin(theta), cos(theta),0], [0,0,1]])
	y = np.dot(DCM,x)
	
	return y

def ecef2ecif(x,t):
	theta = w_earth*t #in radian
	DCM = np.array([[cos(theta), -1*sin(theta), 0], [sin(theta), cos(theta),0],[ 0,0,1]])
	y = np.dot(DCM,x)
	
	return y

def ned2ecef(x,lat,lon):
	v = np.array([-x(2), -x(0), x(1)]) #convert to spherical polar r theta phi
    
    theta = -lat + 90 #in degree, polar angle
    phi = lon #in degree, azimuthal angle
    

    DCM = np.array([[sind(theta)*cosd(phi), cosd(theta)*cosd(phi), -sind(phi)],
    		[sind(theta)*sind(phi), cosd(theta)*sind(phi), cosd(phi)],
    		[cosd(theta), -sind(theta), 0]])
    #for spherical to cartesian

    y = np.dot(DCM,v)
    
    return y