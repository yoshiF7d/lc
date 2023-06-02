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
import AsymmetricGaussianModel
import Baseline
from ParseList import parseList

HEADER_A1 = "LC Chromatogram(Detector A-Ch1)"
HEADER_A2 = "LC Chromatogram(Detector A-Ch2)"
HEADER_B1 = "LC Chromatogram(Detector B-Ch1)"

class ConcParam:
	def __init__(self,conc,param):
		self.conc = conc
		self.param = param
	def __str__(self):
		lst = [self.conc]
		lst.extend(self.param)
		return ','.join([str(l) for l in lst])

class NameConc:
	def __init__(self,name,conc,pos):
		self.name = name
		self.conc = conc
		self.pos = pos

class Chem:
	def __init__(self,name):
		self.name = name
		self.concParams = []
		self.fit = []
		self.pos = 0
	def appendConcParam(self,concParam):
		self.concParams.append(concParam)
	def finalize(self,interpolationOrder):
		self.concParams.sort(key=lambda x:x.conc,reverse=True)
		self.pos = np.average([p.param[0] for p in self.concParams])

		for j in range(4):
			x = [p.conc for p in self.concParams]
			y = [p.param[j] for p in self.concParams]
			if j == 1: 
				# for peak height
				x.append(0)
				y.append(0)
			self.fit.append(np.polyfit(x,y,interpolationOrder))
	def cmax(self):
		return self.concParams[0].conc
	def interpolate(self,conc):
		return np.array(
			[
				np.polyval(self.fit[0],conc),
				np.polyval(self.fit[1],conc),
				np.polyval(self.fit[2],conc),
				np.polyval(self.fit[3],conc),
			]
		)
		
	def print(self):
		print(self.name)
		print(tabulate([[p.conc] + p.param for p in self.concParams],headers=[self.name,'x0','ph','w','sk']))

class Standard:
	def __init__(self,args):
		self.dir = ''
		self.chems = {}
		self.data = []
		self.args = args
	
#	def load(self,file,loader):
#		self.dir = os.path.dirname(file)
#		if os.path.isdir(file):
#			files = os.listdir(file)
#			for f in files:
#				if file.endiwth('.csv'):
#					self.loadstd(os.path.join(file,f),loader)
#		else:
#			self.loadSTD(file)
	
	def load(self,file,loader):
		dir = os.path.dirname(file)
		table = readcsv(file).transpose()
		print(tabulate(table))
		
		if table[-1,0].endswith('.txt'):
			pos = np.full(len(table[0,1:]),None)
		else:
			pos = np.asarray(table[-1,1:],dtype=float)

		for s in table[1:]:
			file = s[0]
			data = loader(os.path.join(dir,file))
			y = data[:,1]
			w = Baseline.calcw(y)
			stdev = np.std(y*w)

			if self.args.noBaselineCorrection:
				data[:,1] -= Baseline.baselineMedian(y)
			else:
				data[:,1] -= Baseline.calcz(y,w)
			
			nameConcs = []
			for i,chemName in enumerate(table[0,1:]):
				conc = float(s[i+1])
				if conc > 0:
					nameConcs.append(NameConc(chemName,conc,pos[i]))
			
			self.setParams(file,data,nameConcs,stdev)
			self.data.append(data)
		
		for chem in self.chems.values():
			chem.finalize(self.args.interpolationOrder)

	def appendChem(self,name,concParam):
		if name not in self.chems:
			self.chems[name] = Chem(name)
		
		self.chems[name].appendConcParam(concParam)
	
	def setParams(self,name,data,nameConcs,stdev):
		x = data[:,0]
		y = data[:,1]
		n = len(data)
		dx = (data[-1,0] - data[0,0]) / n
		w = np.ones(n)
		xi = np.arange(n)

		print(Colors.CYAN + name + Colors.RESET)
		#print(Colors.YELLOW + str(spec) + Colors.RESET)
		peaks,props = find_peaks(y,prominence=self.args.peakProminence*stdev,width=self.args.peakWidth/dx)
		ind = self.peakInclude(x[peaks])
		peaks = peaks[ind]
		rips = props['right_ips'][ind]
		lpis = props['left_ips'][ind]
		#print(peaks)

		if nameConcs[0].pos:		
			ind = []
			for s in nameConcs:
				ind.append(find_nearest(x[peaks],s.pos))
			peaks = peaks[ind]
			width = rips[ind] - lpis[ind]
		else:
			ind = pick(y[peaks],len(nameConcs))
			ind2 = pick(peaks[ind],len(ind))[::-1]
			
			peaks = peaks[ind][ind2]
			width = rips[ind][ind2] - lpis[ind][ind2]

		print(x[peaks])

		pall = []
		blall = []
		brall = []

		for i,p in enumerate(peaks):
			bd = bound(y,p)
			for k in range(n):
				if k < bd[0] or k > bd[1]:
					w[k] = 0
				else:
					w[k] = 1
			
			p0 = [peaks[i],y[peaks[i]],width[i],0]
			p0l = [peaks[i]-width[i],0.5*y[peaks[i]],0.5*width[i],-10]
			p0r = [peaks[i]+width[i],1.5*y[peaks[i]],1.5*width[i],10]
			
			pout = (least_squares(
					AsymmetricGaussianModel.errorAGM,
					x0=p0,
					x_scale='jac',
					#ftol=1e-10,
					jac=AsymmetricGaussianModel.wjacAGM,
					bounds=(p0l,p0r),
					args=(xi,y,w)
				)
			).x

			for po in pout:
				pall.append(po)
			
			blall.extend(p0l)
			brall.extend(p0r)
			pallout = (
				least_squares(self.errorAll,jac=self.jacAll,bounds=(blall,brall),x0=pall,args=(xi,y))
			).x
		
		for i in range(len(peaks)):
			pout = pallout[i*4:(i+1)*4]
			self.appendChem(nameConcs[i].name,ConcParam(nameConcs[i].conc,pout))

	def errorAll(self,p,x,y):
		err = 0
		for i in range(len(p)//4):
			pp = p[4*i:4*(i+1)]
			err += AsymmetricGaussianModel.AGM(x,*pp)
		return y - err
		
	def jacAll(self,p,x,y):
		j =[]
		for i in range(len(p)//4):
			pp = p[4*i:4*(i+1)]
			j.append(AsymmetricGaussianModel.jacAGM(pp,x,y))
		return -np.hstack(j)

	def plotParams(self):
		os.makedirs(self.args.plotParamsDir,exist_ok=True)
		cmax = 0
		for chem in self.chems.values():
			c = chem.cmax()
			if c > cmax:
				cmax = c
			
		for chem in self.chems.values():
			file = os.path.join(self.args.plotParamsDir,chem.name + '.png')
			fig,ax = plt.subplots(2,2)
			label = ['x0','ph','w','sk']

			for j in range(4):
				x = [p.conc for p in chem.concParams]
				y = [p.param[j] for p in chem.concParams]
				x2 = range(0,np.floor(cmax).astype(int),1)
				fit = chem.fit[j]
				
				ax[j//2][j%2].plot(x,y,'bo')
				ax[j//2][j%2].plot(x2,np.polyval(fit,x2),'r-')
				ax[j//2][j%2].set_title(label[j])
				ax[j//2][j%2].set_xlim(0,chem.concParams[0].conc)

				fig.savefig(file)
				plt.close(fig)

	def saveParams(self):
		print(self.args.paramFile)
		with open(self.args.paramFile,'w',newline='') as f:
			s = []
			for chem in self.chems.values():
				row = [','.join([chem.name,'x0','ph','w','sk'])]
				for p in chem.concParams:
					row.append(str(p))
				s.append('\n'.join(row))
			
			f.write('\n'.join(s))
		
	def loadParams(self):
		self.dir = os.path.dirname(self.args.paramFile)
		with open(self.args.paramFile,'r') as f:
			strlist = f.read().split('\n\n')
		
		for s in strlist:
			table = readcsvstr(s)

			for i in range(1,len(table)):
				self.appendChem(table[0][0],ConcParam(float(table[i][0]),np.asarray(table[i][1:],dtype=float)))
			
		for chem in self.chems.values():
			chem.finalize(self.args.interpolationOrder)

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
		
		for i,chem in enumerate(self.chems.values()):
			cmax = chem.concParams[0].conc
			ind = np.floor(chem.pos).astype(int)

			for j in range(10):
				c = cmax*(j+1)/10
				p = chem.interpolate(c)
				ax.plot(x,AsymmetricGaussianModel.AGM(xi,*p),color=clist[i%len(clist)],alpha=0.2)
			if 0 < ind and ind < len(x):
				if self.data:
					ax.text(x[ind],y[ind],chem.name)
				else:
					ax.text(x[ind],p[1],chem.name)
		
		plt.show()
	
	def checkStandards(self):
		fig,ax = plt.subplots()
		for d in self.data:
			x = d[:,0]
			y = d[:,1]
			xi = range(len(x))
			ax.plot(x,y)
		
		y = self.data[0][:,1]
		for chem in self.chems.values():
			ind = np.floor(chem.pos).astype(int)
			ax.text(x[ind],y[ind],chem.name)
		plt.show()

	def print(self):
		for chem in self.chems.values():
			chem.print()
						
	def interpolate(self,chem,conc):
		#return np.array([self.params[ind][0][1](conc),self.params[ind][0][2](conc),self.params[ind][0][3](conc),self.params[ind][0][4](conc)])
		return np.array(
			[
				np.polyval(self.chems[chem].fit[0],conc),
				np.polyval(self.chems[chem].fit[1],conc),
				np.polyval(self.chems[chem].fit[2],conc),
				np.polyval(self.chems[chem].fit[3],conc),
			]
		)

	def eval(self,data,ax):
		chemlen = len(self.chems)
		concs = np.zeros(chemlen+1)
		
		x = data[:,0]
		y = data[:,1]
		xi = np.arange(len(y))
		mask = np.ones(len(y))

		if self.args.peakInclude is not None:
			for lim in self.args.peakInclude:
				mask[np.where(np.logical_and(lim[0] < x, x < lim[1]))] = 0
			mask = 1 - mask

		if self.args.peakExclude is not None:
			for lim in self.args.peakExclude:
				mask[np.where(np.logical_and(lim[0] < x, x < lim[1]))] = 0
		
		#print(mask)

		#initial guess
		for i,chem in enumerate(self.chems.values()):
			c = (data[np.floor(chem.pos).astype(int),1] / chem.concParams[0].param[1]) * chem.concParams[0].conc
			if c > 0:
				concs[i] = c
			else:
				concs[i] = 0

		#add initial guess of x shift to the last of concs
		concs[-1] = 0
		bl = []
		br = []

		for c in concs[:-1]:
			bl.append(0)
			br.append(np.inf)
		
		bl.append(-float(self.args.shiftTolerance))
		br.append(float(self.args.shiftTolerance))

		cout = (least_squares(self.error,x0=concs,x_scale='jac',jac=self.jac,bounds=(bl,br),args=(xi,y*mask))).x
		if ax is not None:
			ax.plot(x,self.model(x,cout),'red',alpha=0.5,linestyle='dashed')
			#ax.plot(data[:,0],self.model(x,1.5*cout),'red',alpha=0.5,linestyle='dashed')
		
		#print("x shift : " + str(cout[-1]))
		return cout[:-1]
		#ax.plot(data[:,0],y,'r',alpha=0.5)
	
	def error(self,concs,x,y):
		return y - self.model(x,concs)
			
	def model(self,x,concs):
		ret = 0
		for i,chem in enumerate(self.chems.values()):
			p = chem.interpolate(concs[i])
			p[0] *= np.exp(concs[-1])
			ret += AsymmetricGaussianModel.AGM(x,*p)
		
		return ret
		
	def jac(self,concs,x,y):
		#dc = 0.01
		res=[]
		las=np.zeros(len(x))
		for i,chem in enumerate(self.chems.values()):
			p = chem.interpolate(concs[i])
			p[0] *= np.exp(concs[-1])
			dpdc = [
				2*concs[i]*chem.fit[0][-3] + chem.fit[0][-2],
				2*concs[i]*chem.fit[1][-3] + chem.fit[1][-2],
				2*concs[i]*chem.fit[2][-3] + chem.fit[2][-2],
				2*concs[i]*chem.fit[3][-3] + chem.fit[3][-2]
			]
			res.append(np.dot(AsymmetricGaussianModel.jacAGM(p,x,y),dpdc))
			las += np.dot(AsymmetricGaussianModel.jacAGM(p,x,y),[1,0,0,0])
		# effect of x shift
		res.append(las)

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

def readcsv(file):
	with open(file) as f:
		data = readcsvstr(f.read())
	return data

def readcsvstr(s):
	row = s.split('\n')
	data = []
	for r in row:
		if r == '':
			continue
		data.append(r.split(','))
	return np.asarray(data,dtype=object)