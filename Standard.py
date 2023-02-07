import argparse
import os
from tabulate import tabulate
import csv
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from scipy.signal import find_peaks
from scipy.signal import peak_widths
from scipy.signal import peak_prominences
from scipy import interpolate
from scipy.optimize import least_squares

from Colors import Colors
from AsymmetricGaussianModel import *
from Baseline import *
from ParseList import parseList

HEADER_A1 = "LC Chromatogram(Detector A-Ch1)"
HEADER_A2 = "LC Chromatogram(Detector A-Ch2)"
HEADER_B1 = "LC Chromatogram(Detector B-Ch1)"

class Standard():
	def __init__(self,args):
		self.dir = ""
		self.data = []
		self.table = []
		self.params = []
		self.peaks = []
		self.args = args
	
	def load(self,stdfile,loader):
		self.dir = os.path.dirname(stdfile)
		with open(stdfile) as f:
			reader = csv.reader(f)
			self.table = [row for row in reader]
		
		temp = [[s[0] for s in self.table]]

		for i,file in enumerate(self.table[0][1:]):
			if file.endswith('.txt'):
				temp.append([s[i+1] for s in self.table])
				self.data.append(loader(os.path.join(self.dir,file)))
			else:
				self.peaks.extend([float(s[i+1]) for s in self.table[1:]])
		
		self.table = np.asarray(temp).transpose()

		if self.args.noBaselineCorrection:
			for data in self.data:
				data[:,1] -= baselineMedian(data[:,1])
		else:
			self.correctStandardBaseline()
		self.setParams()
	
	def setParams(self):
		sane = True
		for s in self.table[1:]:
			self.params.append([[s[0],'x0','ph','w','sk']]) # will be soon replaced by a interpolation object

		for j,s in enumerate(self.table[0][1:]):
			x = self.data[j][:,0]
			y = self.data[j][:,1]
			w = np.ones(len(x))
			xi = range(len(x))
			
			stdev=np.std(y*calcw(y))
			
			print(Colors.CYAN + s + Colors.RESET)
			print('stdev : ' + str(stdev))
			peaks,props = find_peaks(y,prominence=self.args.peakProminence*stdev,width=self.args.peakWidth)
			
			#print(x[peaks])

			ind = self.peakInclude(x[peaks])
			#print(ind)

			peaks = peaks[ind]
			rips = props['right_ips'][ind]
			lpis = props['left_ips'][ind]

			if self.peaks:
				ind = []
				for p in self.peaks:
					ind.append(find_nearest(x[peaks],p))
				peaks = peaks[ind]
				width = rips[ind] - lpis[ind]
			else:
				ind = pick(y[peaks],len(self.params))
				ind2 = pick(peaks[ind],len(ind))[::-1]
				
				peaks = peaks[ind][ind2]
				width = rips[ind][ind2] - lpis[ind][ind2]
			pall = []

			print('peaks : ' + str(x[peaks]))

			for i,p in enumerate(peaks):
				bd = bound(y,p)
				for k in range(len(x)):
					if k < bd[0] or k > bd[1]:
						w[k] = 0
					else:
						w[k] = 1
				
				p0 = [peaks[i],y[peaks[i]],width[i],0]
				p0l = [peaks[i]-width[i],0.5*y[peaks[i]],0.5*width[i],-10]
				p0r = [peaks[i]+width[i],1.5*y[peaks[i]],1.5*width[i],10]
				
				pout = (least_squares(errorAGM,x0=p0,jac=wjacAGM,bounds=(p0l,p0r),args=(xi,y,w))).x
				for p in pout:
					pall.append(p)
			
			if sane:
				pallout = (least_squares(self.errorAll,jac=self.jacAll,x0=pall,args=(xi,y))).x
			else:
				pallout = pall

			if sane:
				for i in range(len(peaks)):
					pout = pallout[i*4:(i+1)*4]
					self.params[i].append([float(self.table[1:][i][j+1]),pout[0],pout[1],pout[2],pout[3]])
				if len(peaks) != len(self.params):
					sane = False
					print(Colors.RED + "warning" + Colors.RESET + " : peak count doesn't match with chemichal count in standard file")
					for i in range(len(peaks),len(self.params)):
						self.params[i].append([0,0,0,1,0])
			else:
				x0_base = np.array([np.average([pp[1] for pp in p[1:]]) for p in self.params])
				x0s = pallout[::4]
				indices=[]
				for i,x in enumerate(x0_base):
					indices.append(np.argmin(np.abs(x0s-x)))
		
				for i in range(len(self.params)):
					k = indices[i]
					pout = pallout[k*4:(k+1)*4]
					self.params[i].append([float(self.table[1:][i][j+1]),pout[0],pout[1],pout[2],pout[3]])

		#for p in self.params:
		#	print(tabulate(p))

		self.makeInterpolation()
	
	def errorAll(self,p,x,y):
		err = 0
		for i in range(len(p)//4):
			pp = p[4*i:4*(i+1)]
			err += AGM(x,*pp)
		return y - err
		
	def jacAll(self,p,x,y):
		j =[]
		for i in range(len(p)//4):
			pp = p[4*i:4*(i+1)]
			j.append(jacAGM(pp,x,y))
		return -np.hstack(j)

	def plotParams(self):
		os.makedirs(self.args.plotParams,exist_ok=True)
		cmax = 0
		for p in self.params:
			if p[1][0] > cmax:
				cmax = p[1][0]
		for p in self.params:
			os.makedirs(os.path.join(self.args.plotParams,p[0][0]),exist_ok=True)
			for j in range(1,len(p[0])):
				x = [p[i][0] for i in range(1,len(p))]
				y = [p[i][j] for i in range(1,len(p))]
				x2 = range(0,np.floor(cmax).astype(int),1)
				fit = p[0][j]
				fig,ax = plt.subplots()
				ax.plot(x,y,'bo')
				ax.plot(x2,fit(x2),'r-')
				ax.set_xlim(0,cmax)
				path = os.path.join(self.args.plotParams,p[0][0],'p'+str(j)+'.png')
				fig.savefig(path)
				plt.close(fig)

	def saveParams(self):
		print(self.args.paramFile)
		with open(self.args.paramFile,'w',newline="") as f:
			writer = csv.writer(f)
			for p in self.params:
				writer.writerow([p[0][0],'x0','ph','w','sk'])
				for row in p[1:]:
					writer.writerow(row)
				f.write('\n')
		
	def loadParams(self):
		self.dir = os.path.dirname(self.args.paramFile)
		with open(self.args.paramFile,'r') as f:
			strlist = f.read().split('\n\n')
		
		for s in strlist:
			param = []
			row = s.split('\n')
			
			if s == "":
				continue
			for r in row:
				param.append(r.split(','))
			for i in range(1,len(param)):
				for j in range(len(param[0])):
					param[i][j] = float(param[i][j])
			
			self.params.append(param)
			
		self.makeInterpolation()

	def checkParams(self):
		fig,ax = plt.subplots()

		if self.data:
			for data in self.data:
				x = data[:,0]
				y = data[:,1]
				ax.plot(x,y,'r',alpha=0.5)

			x = self.data[0][:,0]
			y = self.data[0][:,1]
		else:
			x = np.arange(4801)

		xi = np.arange(len(x))

		if self.args.xlim is not None:
			ax.set_xlim(self.args.xlim)
		if self.args.ylim is not None:
			ax.set_ylim(self.args.ylim)

		clist = ['b','g','c','m','y','k']
		
		for i,pp in enumerate(self.params):
			cmax = pp[1][0]
			ind = np.floor(pp[1][1]).astype(int)
			for j in range(self.args.checkParams):
				c = cmax*(j+1)/self.args.checkParams
				p = self.interpolate(i,c)
				ax.plot(x,AGM(xi,*p),color=clist[i%len(clist)],alpha=0.2)
			if 0 < ind and ind < len(x):
				if self.data:
					ax.text(x[ind],y[ind],pp[0][0])
				else:
					ax.text(x[ind],p[1],pp[0][0])

		plt.show()
	
	def checkStandards(self):
		fig,ax = plt.subplots()
		for d in self.data:
			x = d[:,0]
			y = d[:,1]
			xi = range(len(x))
			ax.plot(x,y)
		
		y = self.data[0][:,1]
		for p in self.params:
			ind = np.floor(p[1][1]).astype(int)
			ax.text(x[ind],y[ind],p[0][0])
		plt.show()

	def print(self):
		for p in self.params:
			print(tabulate(p[1:],headers=[p[0][0],'x0','ph','w','sk']))
	
	def correctStandardBaseline(self):
		data = np.zeros(len(self.data[0]))
		for d in self.data:
			data += d[:,1]
		
		w = calcw(data)
		for k,data in enumerate(self.data):
			data[:,1] -= calcz(data[:,1],w)

	def makeInterpolation(self):
		for p in self.params:
			for j in range(1,len(p[0])):
				x = [p[i][0] for i in range(1,len(p))]
				y = [p[i][j] for i in range(1,len(p))]
				if j == 2:
					x.append(0)
					y.append(0)
				fit = interpolate.interp1d(x,y,kind='linear',fill_value='extrapolate')
				p[0][j] = fit
				
	def interpolate(self,ind,conc):
		return np.array([self.params[ind][0][1](conc),self.params[ind][0][2](conc),self.params[ind][0][3](conc),self.params[ind][0][4](conc)])
	
	def eval(self,data,ax):
		lenp = len(self.params)
		concs = np.zeros(lenp)
		
		y = data[:,1]
		x = range(len(y))
	
		#initial guess
		for i,p in enumerate(self.params):
			x0 = np.average([p[j][1] for j in range(1,len(p))])
			c = (data[np.floor(x0).astype(int),1] / p[1][2]) * p[1][0]
			if c > 0:
				concs[i] = c
			else:
				concs[i] = 0
		
		cout = (least_squares(self.error,x0=concs,jac=self.jac,ftol=1e-10,bounds=(0,np.inf),args=(x,y))).x
		if ax is not None:
			ax.plot(data[:,0],self.model(x,cout),'red',alpha=0.5,linestyle='dashed')
		
		return cout
		#ax.plot(data[:,0],y,'r',alpha=0.5)
	
	def error(self,concs,x,y):
		err = 0
		for i,c in enumerate(concs):
			p = self.interpolate(i,c)
			err += AGM(x,*p)
		
		return y - err
			
	def model(self,x,concs):
		ret = 0
		for i,c in enumerate(concs):
			ret += AGM(x,*self.interpolate(i,c))
		return ret
		
	def jac(self,concs,x,y):
		dc = 0.01
		res=[]
		for i,c in enumerate(concs):
			p = self.interpolate(i,c)
			dpdc = (self.interpolate(i,c + 0.5*dc) - self.interpolate(i,c - 0.5*dc))/dc
			res.append(np.dot(jacAGM(p,x,y),dpdc))
		return - np.transpose(res)

	def peakInclude(self,lst):
		peak_set = set()
		for i,x in enumerate(lst):
			if self.args.peakInclude is not None:
				isin = False
				for lim in self.args.peakInclude:
					if lim[0] < x and x < lim[1]:
						isin = True
						break
			else:
				isin = True
			if isin and self.args.peakExclude is not None:
				for lim in self.args.peakExclude:
					if lim[0] < x and x < lim[1]:
						isin = False
						break
			if isin:
				peak_set.add(i)
		
		return np.array(sorted(list(peak_set)))

def bound(ylist, peak):
	xmin = 0
	xmax = len(ylist)
	left_bound = right_bound = peak
	dy = 0
	left_bound -= 1
	right_bound += 1
	while left_bound >= xmin:
		dy = 0
		i = 0
		while dy >= 0 and i < 5:
			if left_bound - i - 1 >= xmin:
				dy += ylist[left_bound - i] - ylist[left_bound - i - 1]
				i += 1
			else:
				dy = -1
				break
		if dy < 0:
			break
		left_bound -= 1

	while right_bound < xmax:
		dy = 0
		i = 0
		while dy >= 0 and i < 5:
			if right_bound + i + 1 < xmax:
				dy += ylist[right_bound + i] - ylist[right_bound + i + 1]
				i += 1
			else:
				dy = -1
				break
		if dy < 0:
			break
		right_bound += 1

	return [left_bound, right_bound]

def pick(data,count):
	return np.argsort(data)[::-1][:count]

def find_nearest(lst,val):
	lst = np.asarray(lst)
	return np.abs(lst-val).argmin()