#!/usr/local/bin/python
from math import *
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import os
import errno
import numpy as np

## TimeIntegrationSolvers.py
## Remi Carmigniani
#
# Example : y'(t) = -cos(t)*y(t) with IC y(t) = 1 
# Solution is y(t) = y(0)*Exp(-sin(t))
#More generally y'(t) = f(y,t)
#Reminder on the schemes : 
#	Euler Explicit : y(t+dt) = y(t)+dt*f(y(t),t)
# 	2nd Runge Kutta : y(t+dt) = 

directory = 'result'
if not os.path.exists(directory):
    os.makedirs(directory)


##################################################################################################################
############################			Useful functions		      ############################
##################################################################################################################
def f(y,t):
	return -cos(t)*y

def euler(y,t,dt):
	return y+dt*f(y,t)

def RK2(y,t,dt):
	return y+dt*f(y+0.5*dt*f(y,t),t+0.5*dt)
	

def RK4(y,t,dt):
	k1=f(y,t)
	k2 = f(y+.5*dt*k1,t+.5*dt)
	k3 = f(y+.5*dt*k2,t+.5*dt)
	k4 = f(y+dt*k3,t+dt)
	return y+dt*(1./6.*k1+1./3.*k2+1./3.*k3+1./6.*k4)


def exact(t):
	return u0*exp(-sin(t))
##################################################################################################################
############################				End			      ############################
##################################################################################################################



## Discretization parameter 
dtTable = [0.5,0.2,0.1,0.01,0.001,0.0001,0.00001]
tend = 20.
u0 = 1.
errorEuler = []
error2RK = []
error4RK = []
for i in range(len(dtTable)):
	t=0
	dt = dtTable[i]
	solEuler = []
	sol2RK = []
	sol4RK = []
	time = []
	solEuler.append(u0)
	sol2RK.append(u0)
	sol4RK.append(u0)
	time.append(t)
	while t< tend:
		solEuler.append(euler(solEuler[-1],t,dt))
		sol2RK.append(RK2(sol2RK[-1],t,dt))
		sol4RK.append(RK4(sol4RK[-1],t,dt))
		t=t+dt
		time.append(t)


	plt.plot(time,solEuler,label='Euler')
	plt.plot(time,sol2RK,label='2nd RK')
	plt.plot(time,sol4RK,label='4th RK')
	plt.plot(time, [exact(time[ii]) for ii in range(len(time))],label = 'exact',ls='--')
	plt.legend(loc=1)
	plt.title('Solution for dt='+repr(dt))
	plt.savefig('Solutiondt_'+repr(i)+'.png')
	plt.close()
	
	errorEuler.append(sqrt((solEuler[-1]-exact(time[-1]))**2))
	error2RK.append(sqrt((sol2RK[-1]-exact(time[-1]))**2))
	error4RK.append(sqrt((sol4RK[-1]-exact(time[-1]))**2))


plt.loglog(dtTable, errorEuler,label='Euler')
plt.loglog(dtTable, error2RK,label='2nd RK')
plt.loglog(dtTable, error4RK,label='4th RK')
plt.legend(loc=1)
plt.title('Error')
plt.savefig('Error.png')





print 'Simulation Completed without error'


