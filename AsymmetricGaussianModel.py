import numpy as np
from scipy.special import erf

#https://stackoverflow.com/questions/62730830/asymmetric-gaussian-fit-in-python
def AGM(x,x0,ph,w,sk):
	e = (x-x0)/w
	return ph*(1 + erf(sk*e))*np.exp(-e*e)

def errorAGM(p,x,y,w):
	return w*(AGM(x,*p) - y)
	
def wjacAGM(p,x,y,w):
	return np.transpose(w*np.transpose(jacAGM(p,x,y)))
	
def jacAGM(p,x,y):
	e = (x-p[0])/p[2]
	f = 1+erf(e*p[3])
	g = np.exp(-e*e*p[3]*p[3])/np.sqrt(np.pi)
	h = p[1]*(-p[3]*g + e*f)/p[2]
	j0 = 2*h
	j1 = f
	j2 = 2*e*h
	j3 = 2*e*g*p[1]

	return np.transpose([j0,j1,j2,j3]*np.exp(-e*e))