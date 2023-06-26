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

class Chem:
	def __init__(self,name,conc,param):
		self.name = name
		self.concs = [conc]
		self.params = [param]
		self.fit = []
		self.pos = param[0]
	
	def append(self,conc,param):
		self.concs.append(conc)
		self.params.append(param)
		self.pos = 0.5 * (self.pos + param[0])

	def finalize(self):
		self.concs = np.asarray(self.concs)
		self.params = np.asarray(self.params)
		inds = (-self.concs).argsort()
		self.concs = self.concs[inds]
		self.params = self.params[inds]
		self.fit = []
		self.fit.append(np.polyfit(self.concs,self.params[:,0],2))
		self.fit.append(np.polyfit(np.append([0],self.concs),np.append([0],self.params[:,1]),2))
		self.fit.append(np.polyfit(self.concs,self.params[:,2],2))
		self.fit.append(np.polyfit(self.concs,self.params[:,3],2))		

	def interpolate(self,conc):
		return np.array(
			[
				np.polyval(self.fit[0],conc),
				np.polyval(self.fit[1],conc),
				np.polyval(self.fit[2],conc),
				np.polyval(self.fit[3],conc),
			]
		)

		# p(c) = a[0]*c*c + a[1]*c + a[2]
		# dpdc = 2*c*a[0] + a[1] 

	def dpdc(self,conc):
		return [
			2*conc*self.fit[0][0] + self.fit[0][1],
			2*conc*self.fit[1][0] + self.fit[1][1],
			2*conc*self.fit[2][0] + self.fit[2][1],
			2*conc*self.fit[3][0] + self.fit[3][1]
		]

	def print(self):
		print(self.name)
		print(tabulate([[c] + p for c,p in zip(self.concs,self.params)]))

class Standard:
	def __init__(self,args):
		self.dir = ''
		self.chems = {}
		self.data = []
		self.args = args
	
	def appendChem(self,name,conc,param):
		if name not in self.chems:
			self.chems[name] = Chem(name,conc,param)
		else:
			self.chems[name].append(conc,param)
	
	def getpos(self,names):
		pos = []
		for name in names:
			if name in self.chems:
				pos.append(self.chems[name].pos)
			else:
				pos.append(None)
		return np.asarray(pos)
			
	def load(self,file,loader):
		dir = os.path.dirname(file)
		table = readcsv(file).transpose()
		print(tabulate(table))
		
		names = np.asarray(table[0,1:])
		concs = np.zeros(len(names))
		
		if not table[-1,0].endswith('.txt'):
			pos = np.asarray(table[-1,1:],dtype=float)
			lst = table[1:-1]
		else:
			pos = self.getpos(names)
			lst = table[1:]

		for s in lst:
			file = s[0]
			data = loader(os.path.join(dir,file))
			y = data[:,1]
			w = Baseline.calcw(y)
			stdev = np.std(y*w)

			if self.args.noBaselineCorrection:
				data[:,1] -= Baseline.baselineMedian(y)
			else:
				data[:,1] -= Baseline.calcz(y,w)
			
			for i,name in enumerate(names):
				concs[i] = tofloat(s[i+1])

			ind = concs > 0

			Names = names[ind]
			Concs = concs[ind]
			Pos = pos[ind]

			self.setParams(file,data,Names,Concs,Pos,stdev)
			
			self.data.append(data)
		
		for chem in self.chems.values():
			chem.finalize()
	
	def setParams(self,file,data,names,concs,pos,stdev):
		x = data[:,0]
		y = data[:,1]
		n = len(data)
		dx = (data[-1,0] - data[0,0]) / n
		w = np.ones(n)
		xi = np.arange(n)

		peaks,props = find_peaks(y,prominence=self.args.peakProminence*stdev,width=self.args.peakWidth/dx)
		inds = self.peakInclude(x[peaks])
		peaks,rips,lips = peaks[inds],props['right_ips'][inds],props['left_ips'][inds]
		inds = np.arange(len(peaks))
		
		Inds = []
		Peaks = []
		Widths = []

		for k,p in enumerate(pos):
			if p:
				#print(names[k])
				i = find_nearest(x[peaks],p)
				Inds.append(i)
				Peaks.append(peaks[i])
				Widths.append(rips[i] - lips[i]) 
			else:
				Peaks.append(None)
				Widths.append(None)
		
		cp = len(Inds) #count of peaks with known pos 
		cuk = len(names) - cp #count of peaks with unknown pos

		ind = complement(inds,Inds)
		peaks,rips,lips = peaks[ind],rips[ind],lips[ind]

		ind = pick(y[peaks],cuk)
		ind2 = pick(peaks[ind],len(ind))[::-1]
		
		peaks = peaks[ind][ind2]
		width = rips[ind][ind2] - lips[ind][ind2]

		j = 0
		for i in range(len(names)):
			if Peaks[i] is None:
				Peaks[i] = peaks[j]
				Widths[i] = width[j]
				j = j + 1

		pall = []
		blall = []
		brall = []

		print(tabulate([[file] +[f'stdev : {stdev:.1g}'] + ['peaks :'] + [f'({x[p]:.1f} {y[p]:.2g})' for p in Peaks]],tablefmt="plain",floatfmt=".1f"))

		for peak,width in zip(Peaks,Widths):
			bd = bound(y,peak)
			for k in range(n):
				if k < bd[0] or k > bd[1]:
					w[k] = 0
				else:
					w[k] = 1
			
			p0 = [peak,y[peak],width,0]
			p0l = [peak-width,0.5*y[peak],0.5*width,-10]
			p0r = [peak+width,1.5*y[peak],1.5*width,10]
			
			pout = (
				least_squares(
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
		
		for i in range(len(Peaks)):
			pout = pallout[i*4:(i+1)*4]
			self.appendChem(names[i],concs[i],pout)

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
		os.makedirs(self.args.plotDir,exist_ok=True)
		cmax = 0
		for chem in self.chems.values():
			c = chem.concs[0]
			if c > cmax:
				cmax = c
			
		for chem in self.chems.values():
			file = os.path.join(self.args.plotDir,chem.name + '.png')
			fig,ax = plt.subplots(2,2)
			label = ['x0','ph','w','sk']

			x = chem.concs
			x2 = np.linspace(0,np.floor(cmax).astype(int),100)

			for j in range(4):
				y = chem.params[:,j]
				fit = chem.fit[j]
				ax[j//2][j%2].plot(x,y,'bo')
				ax[j//2][j%2].plot(x2,np.polyval(fit,x2),'r-')
				ax[j//2][j%2].set_title(label[j])
				ax[j//2][j%2].set_xlim(0,chem.concs[0])

				fig.savefig(file)
				plt.close(fig)
			
		fig,ax = plt.subplots()
		self.checkParams(ax)
		fig.savefig(os.path.join(self.args.plotDir,'fit.png'))

	def saveParams(self):
		print(self.args.paramFile)
		with open(self.args.paramFile,'w',newline='') as f:
			s = []
			for chem in self.chems.values():
				row = [','.join([chem.name,'x0','ph','w','sk'])]
				for conc,param in zip(chem.concs,chem.params):
					row.append(','.join([str(l) for l in np.append([conc],param)]))
				s.append('\n'.join(row))
			
			f.write('\n\n'.join(s))
		
	def loadParams(self):
		self.dir = os.path.dirname(self.args.paramFile)
		with open(self.args.paramFile,'r') as f:
			strlist = f.read().split('\n\n')
		
		for s in strlist:
			table = readcsvstr(s)
			for i in range(1,len(table)):
				self.appendChem(table[0][0],float(table[i][0]),np.asarray(table[i][1:],dtype=float))

		for chem in self.chems.values():
			for i,p in enumerate(chem.params):
				chem.params[i][0] *= self.args.shift
			chem.finalize()

	def checkParams(self,ax=None):
		show = False
		if not ax:
			show = True
			fig,ax = plt.subplots()

		if self.data:
			for data in self.data:
				x = data[:,0]
				y = data[:,1]
				ax.plot(x,y,'r',alpha=0.5)

			x = self.data[0][:,0]
		else:
			x = np.arange(4801)

		xi = np.arange(len(x))

		if self.args.xlim is not None:
			ax.set_xlim(self.args.xlim)
		if self.args.ylim is not None:
			ax.set_ylim(self.args.ylim)

		clist = ['b','g','c','m','y','k']
		
		for i,chem in enumerate(self.chems.values()):
			cmax = chem.concs[0]
			ind = np.floor(chem.pos).astype(int)

			for j in range(10):
				c = cmax*(j+1)/10
				p = chem.interpolate(c)
				ax.plot(x,AsymmetricGaussianModel.AGM(xi,*p),color=clist[i%len(clist)],alpha=0.2)
			if 0 < ind and ind < len(x):
				ax.text(x[ind],chem.params[0][1],chem.name)
		
		if show:
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

	def eval(self,data,ax):
		chemlen = len(self.chems)
		concs = np.zeros(chemlen)
		
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

		#initial guess
		for i,chem in enumerate(self.chems.values()):
			c = (data[np.floor(chem.pos).astype(int),1] / chem.params[0][1]) * chem.concs[0]
			if c > 0:
				concs[i] = c
			else:
				concs[i] = 0

		bl = []
		br = []

		for c in concs:
			bl.append(0)
			br.append(np.inf)
		
		if self.args.shiftTolerance is not None and self.args.shift == 1:
			concs = np.append(concs,[0])
			bl.append(-float(self.args.shiftTolerance))
			br.append(float(self.args.shiftTolerance))
			cout = (least_squares(self.errorShift,x0=concs,x_scale='jac',jac=self.jacShift,bounds=(bl,br),args=(xi,y*mask))).x
			self.args.shift = np.exp(cout[-1])
			self.applyShift()
			bl = bl[:-1]
			br = br[:-1]
			concs = concs[:-1]

		cout = (least_squares(self.error,x0=concs,x_scale='jac',jac=self.jac,bounds=(bl,br),args=(xi,y*mask))).x
		if ax is not None:
			ax.plot(x,self.model(xi,cout),'red',alpha=0.5,linestyle='dashed')
		
		return cout
	
	def error(self,concs,x,y):
		return y - self.model(x,concs)
	
	def errorShift(self,concs,x,y):
		return y - self.modelShift(x,concs)

			
	def model(self,x,concs):
		ret = 0
		for i,chem in enumerate(self.chems.values()):
			p = chem.interpolate(concs[i])
			ret += AsymmetricGaussianModel.AGM(x,*p)
		
		return ret

	def modelShift(self,x,concs):
		ret = 0
		for i,chem in enumerate(self.chems.values()):
			p = chem.interpolate(concs[i])
			p[0] *= np.exp(concs[-1])
			ret += AsymmetricGaussianModel.AGM(x,*p)
		
		return ret
	
	def jac(self,concs,x,y):
		#dc = 0.01
		res=[]
		for i,chem in enumerate(self.chems.values()):
			p = chem.interpolate(concs[i])
			dpdc = chem.dpdc(concs[i])
			res.append(np.dot(AsymmetricGaussianModel.jacAGM(p,x,y),dpdc))

		return - np.transpose(res)
		
	def jacShift(self,concs,x,y):
		#dc = 0.01
		res=[]
		las=np.zeros(len(x))
		for i,chem in enumerate(self.chems.values()):
			p = chem.interpolate(concs[i])
			p[0] *= np.exp(concs[-1])
			dpdc = chem.dpdc(concs[i])
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

	def applyShift(self):
		print('set x shift multiplier : ' + str(self.args.shift))
		for chem in self.chems.values():
			for i,p in enumerate(chem.params):
				chem.params[i,0] = p[0]*self.args.shift
			chem.finalize()


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

def complement(e_all, *args):
    e_all = set(e_all)
    for e in args:
        e_all -= set(e)
    return np.asarray(list(e_all))

def tofloat(s):
	try:
		f = float(s)
	except:
		f = 0
	return f