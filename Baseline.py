import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.special import expit

def baselineMedian(y):
	return np.full(len(y),np.median(y))	

#https://stackoverflow.com/questions/29156532/python-baseline-correction-library
def baseline(data):
	w = calcw(data)
	return calcz(data,w)

def calcw(data):
	dif = np.diff(data,append=0)
	std = np.std(dif)
	hist,x = np.histogram(dif,bins=100,range=(-std,std),density=True)
	s = 1.0/(np.max(hist)*np.sqrt(np.pi))
	
	w = np.exp(-dif*dif/(s*s))
	z = np.array(range(101))
	z = np.exp(-(z-50)*(z-50)/(5*5))
	w = np.convolve(w,z,mode='same')
	#fig,ax = plt.subplots()
	x = np.convolve(x,[0.5,0.5],mode='valid')
	#ax.plot(range(len(w)),w*1e+4)
	#ax.plot(x,max(hist)*np.exp(-x*x/(s*s)))
	#ax.plot(range(len(data)),data)
	#ax.add_patch(patches.Rectangle((0,-0.5*s),len(dif),s,facecolor='orange',alpha=0.3))
	#ax.set_ylim([-500,500])
	#plt.show()
	#exit()
	return w/w.max()

def calcz(data,w):
	y = data
	l = len(y)
	D = sparse.diags([1,-2,1],[0,-1,-2], shape=(l,l-2))
	D = 1e+7*D.dot(D.transpose())
	W = sparse.spdiags(w,0,l,l)
	for i in range(10):
		W.setdiag(w)
		Z = W + D
		z = spsolve(Z, w*y)
		d = y-z
		dn = d[d<0]
		m = np.mean(dn)
		s = np.std(dn)
		wnew = w*expit(2 * ((2*s - m)-d)/s)
		w = wnew
	return z