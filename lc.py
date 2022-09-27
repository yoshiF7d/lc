import argparse
import os
from tabulate import tabulate
import csv
import numpy as np
import re
import matplotlib.pyplot as plt
#import matplotlib.patches as patches

from scipy.signal import find_peaks
from scipy.signal import peak_widths
from scipy.signal import peak_prominences
from scipy import interpolate
from scipy.optimize import least_squares
from scipy.special import erf
from scipy.special import expit
from scipy import sparse
from scipy.sparse.linalg import spsolve

os.system('')

HEADER = "LC Chromatogram(Detector B-Ch1)"
#HEADER = "LC Chromatogram(Detector A-Ch2)"
LAMBDA = 105

class colors():
	RESET = '\033[0m'
	GRAY = '\033[90m'
	RED = '\033[91m'
	GREEN = '\033[92m'
	YELLOW = '\033[93m'
	BLUE = '\033[94m'
	MAGENTA = '\033[95m'
	CYAN = '\033[96m'
	DARKRED = '\033[31m'
	DARKGREEN = '\033[32m'
	DARKBLUE = '\033[34m'
	MARINE = '\033[35m'
	DARKCYAN = '\033[36m'

class LCData():
	def __init__(self,file):
		self.content = []
		self.load(file)

	def load(self,file):
		with open(file,encoding='Shift_JIS') as f:
			strlist = f.read().split('\n\n')

		for s in strlist:
			t = s.strip()
			if t:
				self.content.append(LCSubData(t))

	def print(self):
		for c in self.content:
			c.print()
			print("\n")

	def query(self,key):
		for c in self.content:
			if key in c.name:
				if c.type == "LC Chromatogram" or c.type == "LC Status Trace":
					c.readT(c.rest,False)
				break

		return c.data
	
class LCSubData():
	def __init__(self,strin):
		self.name = None
		self.table = None
		self.data = None
		self.type = None
		self.load(strin,headerOnly=True)

	def load(self,strin,headerOnly=False):
		self.name = strin.split('\n',maxsplit=1)
		self.rest = None
		if len(self.name) > 1:
			self.rest = self.name[1]
		
		self.name = self.name[0].strip('[]')
		
		if not self.rest:
			return

		if bool(re.match("Peak Table",self.name)):
			self.type = "Peak Table"
			self.readR(self.rest)
		elif bool(re.match("Compound Results",self.name)):
			self.type = "Compound Results"
			self.readR(self.rest)
		elif bool(re.match("Group Results",self.name)):
			self.type = "Group Results"
			self.readR(self.rest)
		elif bool(re.match("LC Chromatogram",self.name)):
			self.type = "LC Chromatogram"
			self.readT(self.rest,headerOnly)
		elif bool(re.match("LC Status Trace",self.name)):
			self.type = "LC Status Trace"
			self.readT(self.rest,headerOnly)
		else:
			self.read(self.rest)

		if len(self.table) > 1:
			if type(self.table[1]) is list:
				return
			
		self.table = [self.table]
	
	def print(self):
		#print(colors.RED + str(self.type) + colors.RESET)
		print(colors.GREEN + self.name + colors.RESET)
		#print(colors.RED + str(self.table) + colors.RESET)
		print(colors.YELLOW + tabulate(self.table) + colors.RESET)
		
		if	self.data is not None:
			if type(self.data) is np.ndarray:
				print(self.data)
			else:
				print(tabulate(self.data))
	
	def read(self,strin):
		lines = strin.splitlines()
		for j,line in enumerate(lines):
			lines[j] = line.split("\t")
		
		if len(lines) == 1:
			self.table = [lines]
		else:
			self.table = lines

	def readR(self,strin):
		lines = strin.splitlines()
		for j,line in enumerate(lines):
			lines[j] = line.split("\t")
		
		if len(lines) == 1:
			self.table = [lines[0]]
		else:
			self.table = lines[0]

		self.data = lines[1:] 
		
	def readT(self,strin,headerOnly=True):
		self.table = []
		s = strin	
	
		while s:
			rest = None
			lines = s.split("\n",maxsplit=1)
			if len(lines) > 1:
				rest = lines[1]
			
			if "R.Time" in lines[0]:
				if not headerOnly:
					 self.loadData(lines[1])
				break
				
			self.table.append(lines[0].split("\t"))
			s = rest	

		#self.data = np.loadtxt(lines[1])[:5]

	def loadData(self,strin):
		self.data = []
		lines = strin.split("\n")
		for s in lines:
			column = s.split("\t")
			cdata = []
			for c in column:
				cdata.append(float(c))
			self.data.append(cdata)
		
		self.data = np.array(self.data)

class STDT():
	def __init__(self):
		self.dir = ""
		self.data = []
		self.table = []
		self.params = []
	
	def load(self,stdfile,noBLE=False):
		self.dir = os.path.dirname(stdfile)
		with open(stdfile) as f:
			reader = csv.reader(f)
			self.table = [row for row in reader]

		for file in list(map(lambda s:os.path.join(self.dir,s),self.table[0][1:])):
			if file.endswith('.txt'):
				#data = LCData(file).query(args.header)
				#data[:,1] -= baseline(data[:,1])
				self.data.append(LCData(file).query(args.header))
		
		if noBLE:
			for d in data:
				d[:,1] -= baselineMedian(data[:,1])
		else:
			self.correctSTDBL()
		self.setParams()
	
	def setParams(self):
		for s in self.table[1:]:
			self.params.append([[s[0],'x0','ph','w','sk']]) # will be soon replaced by a interpolation object

		for j,s in enumerate(self.table[0][1:]):
			x = self.data[j][:,0]
			y = self.data[j][:,1]
			w = np.ones(len(x))
			xi = range(len(x))
			
			stdev=np.std(y)
			peaks,props = find_peaks(y,prominence=0.1*stdev,width=15)
			ind = pick(y[peaks],len(self.table)-1)
			ind2 = pick(peaks[ind],len(ind))[::-1]
			
			peaks = peaks[ind][ind2]
			width = props["right_ips"][ind][ind2] - props["left_ips"][ind][ind2]
			pall = []

			for i,p in enumerate(peaks):
				bd = bound(y,p)
				#k0 = (bd[0] + bd[1])*0.5
				#w0 = (bd[1] - bd[0])*0.5
				for k in range(len(x)):
					#w[k] = 0.5*(1-np.tanh(WP0*(np.abs(k - k0)-WP1*w0)))
					if k < bd[0] or k > bd[1]:
						w[k] = 0
					else:
						w[k] = 1
				
				p0 = [peaks[i],y[peaks[i]],width[i],0]
				pout = (least_squares(errorAGM,x0=p0,jac=wjacAGM,args=(xi,y,w))).x
				for p in pout:
					pall.append(p)
				#self.params[i].append([float(self.table[1:][i][j+1]),pout[0],pout[1],pout[2],pout[3]])
			
			pallout = (least_squares(self.errorWhole,jac=self.jacEW,x0=pall,args=(xi,y))).x
			
			for i,p in enumerate(peaks):
				pout = pallout[i*4:(i+1)*4]
				self.params[i].append([float(self.table[1:][i][j+1]),pout[0],pout[1],pout[2],pout[3]])
		
		self.makeintp()
	
	def errorWhole(self,p,x,y):
		err = 0
		for i in range(len(self.params)):
			pp = p[4*i:4*(i+1)]
			err += asymGaussModel(x,*pp)
		return y - err
		
	def jacEW(self,p,x,y):
		j =[]
		for i in range(len(self.params)):
			pp = p[4*i:4*(i+1)]
			j.append(jacAGM(pp,x,y))
		return -np.hstack(j)

	def autoFit(self):
		pass

	def plotParams(self,pdir):
		os.makedirs(pdir,exist_ok=True)
		cmax = 0
		for p in self.params:
			if p[1][0] > cmax:
				cmax = p[1][0]
		for p in self.params:
			os.makedirs(os.path.join(pdir,p[0][0]),exist_ok=True)
			for j in range(1,len(p[0])):
				x = [p[i][0] for i in range(1,len(p))]
				y = [p[i][j] for i in range(1,len(p))]
				x2 = range(0,cmax,1)
				fit = p[0][j]
				fig,ax = plt.subplots()
				ax.plot(x,y,'bo')
				ax.plot(x2,fit(x2),'r-')
				ax.set_xlim(0,cmax)
				path = os.path.join(pdir,p[0][0],'p'+str(j)+'.png')
				fig.savefig(path)
				plt.close(fig)

	def saveParams(self,pfile):
		with open(pfile,'w',newline="") as f:
			writer = csv.writer(f)
			for p in self.params:
				writer.writerow([p[0][0],'x0','ph','w','sk'])
				for row in p[1:]:
					writer.writerow(row)
				f.write('\n')
		
	def loadParams(self,pfile):
		self.dir = os.path.dirname(pfile)
		with open(pfile,'r') as f:
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
			
		self.makeintp()

	def checkParams(self,count,xlim,ylim):
		fig,ax = plt.subplots()
		for data in self.data:
			x = data[:,0]
			y = data[:,1]
			xi = range(len(x))
			ax.plot(x,y,'r',alpha=0.5)
		
		if xlim is not None:
			ax.set_xlim(xlim)
		if ylim is not None:
			ax.set_ylim(ylim)

		clist = ['b','g','c','m','y','k']
		
		for i,pp in enumerate(self.params):
			cmax = pp[1][0]
			ind = np.floor(pp[1][1]).astype(int)
			ax.text(x[ind],self.data[0][ind,1],pp[0][0])
			for j in range(count):
				c = cmax*(j+1)/count
				p = stdt.intp(i,c)
				ax.plot(x,asymGaussModel(xi,*p),color=clist[i%len(clist)],alpha=0.2)
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
	
	def correctSTDBL(self):
		data = np.zeros(len(self.data[0]))
		for d in self.data:
			data += d[:,1]
		
		w = calcw(data)
		#fig,ax = plt.subplots(len(self.data))
		#for d in self.data:
		#fig,ax = plt.subplots(2)
		#for data in self.data:
		#	ax[0].plot(range(len(data)),data[:,1])
		#ax[0].set_ylim([-10000,10000])
		
		for k,data in enumerate(self.data):
			data[:,1] -= calcz(data[:,1],w)
			#ax[k].plot(range(len(y)),z,color='b',alpha=0.5)
			#ax[k].plot(range(len(y)),y,color='r',alpha=0.5)
			#ax[k].set_ylim([-10000,10000])
		#for data in self.data:
		#	ax[1].plot(range(len(data)),data[:,1])
		#ax[1].set_ylim([-10000,10000])
		#plt.show()
		#exit()
		
	def makeintp(self):
		for p in self.params:
			for j in range(1,len(p[0])):
				x = [p[i][0] for i in range(1,len(p))]
				y = [p[i][j] for i in range(1,len(p))]
				if j == 2:
					x.append(0)
					y.append(0)
				fit = interpolate.interp1d(x,y,kind='linear',fill_value='extrapolate')
				p[0][j] = fit
				
	def intp(self,ind,conc):
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
			ax.plot(data[:,0],self.model(x,cout),'b',alpha=0.5,linestyle='dashed')
		
		return cout
		#ax.plot(data[:,0],y,'r',alpha=0.5)
	
	def error(self,concs,x,y):
		err = 0
		for i,c in enumerate(concs):
			p = self.intp(i,c)
			err += asymGaussModel(x,*p)
		
		return y - err
			
	def model(self,x,concs):
		ret = 0
		for i,c in enumerate(concs):
			ret += asymGaussModel(x,*self.intp(i,c))
		return ret
		
	def jac(self,concs,x,y):
		dc = 0.01
		res=[]
		for i,c in enumerate(concs):
			p = self.intp(i,c)
			dpdc = (self.intp(i,c + 0.5*dc) - self.intp(i,c - 0.5*dc))/dc
			res.append(np.dot(jacAGM(p,x,y),dpdc))
		return - np.transpose(res)
		
def bound(ylist,peak):
	lb = rb = peak
	xmin = 0
	xmax = len(ylist)
	dy = 0
	lb -= 1
	rb += 1
	while lb >= xmin:
		dy = 0
		for i in range(5):
			if lb-i-1 >= xmin:
				dy += ylist[lb-i] - ylist[lb-i-1]
			else:
				dy = -1
				break
		
		if dy < 0:
			break
			
		lb -= 1

	while rb < xmax:
		dy = 0
		for i in range(5):
			if rb+i+1 < xmax:
				dy += ylist[rb+i] - ylist[rb+i+1]
			else:
				dy = -1
				break
		
		if dy < 0:
			break
			
		rb += 1

	return [lb,rb]
	
def pick(data,count):
	return np.argsort(data)[::-1][:count]

#https://stackoverflow.com/questions/62730830/asymmetric-gaussian-fit-in-python
def asymGaussModel(x,x0,ph,w,sk):
	e = (x-x0)/w
	return ph*(1 + erf(sk*e))*np.exp(-e*e)

def errorAGM(p,x,y,w):
	return w*(asymGaussModel(x,*p) - y)
	
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

#https://stackoverflow.com/questions/4623446/how-do-you-sort-files-numerically
def tryint(s):
	try:
		return int(s)
	except:
		return s

def alnum_key(s):
	return [tryint(c) for c in re.split('([0-9]+)',s)]

def parselim(limstr):
	l0,l1 = limstr.split(",")
	l0 = l0.lstrip('([').lstrip()
	l1 = l1.rstrip(')]').rstrip()
	return [float(l0),float(l1)]

parser = argparse.ArgumentParser(description='LC processor')
parser.add_argument('filein',nargs='?',help='Input file or directory')
parser.add_argument('fileout',nargs='?',help='Output file')
parser.add_argument('-s','--stdFile',help="Standard file")
parser.add_argument('-p','--paramFile',help="Parameter file")
parser.add_argument('--plotDir',help="Plot Directory")
parser.add_argument('--plotParamsDir',help='Plot parameters Directory')
parser.add_argument('--checkParams',action='store_const',const=10,help='check parameters (Standard file is required)')
parser.add_argument('--checkStandards',action='store_true',help='check standards data (Standard file is required)')
parser.add_argument('--header',default=HEADER,help='header name to query in a datafile')
parser.add_argument('--noBLE',default=False,action='store_true',help='no baseline estimation')
parser.add_argument('--xlim',help='x region of plot')
parser.add_argument('--ylim',help='y region of plot')
parser.add_argument('--label',help='show chemical species label',action='store_true')

stdt = STDT()
args = parser.parse_args()
if args.xlim is not None:
	args.xlim = parselim(args.xlim)
if args.ylim is not None:
	args.ylim = parselim(args.ylim)

if args.stdFile is not None:
	stdt.load(args.stdFile,args.noBLE)
	print(colors.YELLOW + "Standard Data Table" + colors.RESET)
	print(tabulate(stdt.table))
	if args.paramFile is not None:
		stdt.saveParams(args.paramFile)
	if args.checkParams:
		stdt.checkParams(int(args.checkParams),args.xlim,args.ylim)
		exit()
	if args.checkStandards:
		stdt.checkStandards()
		exit()
	if args.plotParamsDir is not None:
		stdt.plotParams(args.plotParamsDir)
		exit()
elif args.paramFile is not None:
	stdt.loadParams(args.paramFile)
	if args.plotParamsDir is not None:
		stdt.plotParams(args.plotParamsDir)
		exit()
else:
	lcdata = LCData(args.filein)
	data = lcdata.query(args.header)
	if args.noBLE:
		data[:,1] -= baselineMedian(data[:,1])
	else:
		data[:,1] -= baseline(data[:,1])
	
	x = data[:,0]
	y = data[:,1]

	fig,ax = plt.subplots()
	ax.plot(x,y,'r')
	if args.xlim is not None:
		ax.set_xlim(args.xlim)
	if args.ylim is not None:
		ax.set_ylim(args.ylim)
	plt.show()
	#print('stdFile or paramFile is required')
	exit()	

if args.filein is not None:
	if os.path.isdir(args.filein):
		#multiple file mode	
		files = os.listdir(args.filein)
		files.sort(key=alnum_key)
		row = ['File Name']+[p[0][0] for p in stdt.params]
		fmtstr = "{:<20}"
		for p in stdt.params:
			fmtstr += " {:<15}"
		print(fmtstr.format(*row))
		output = []
		output.append(row)
		if args.plotDir is not None:
			os.makedirs(args.plotDir,exist_ok=True)
		for file in files:
			if file.endswith('.txt'):
				lcdata = LCData(os.path.join(args.filein,file))
				data = lcdata.query(args.header)
				if args.noBLE:
					data[:,1] -= baselineMedian(data[:,1])
				else:
					data[:,1] -= baseline(data[:,1])
			
				if args.plotDir is not None:
					fig,ax = plt.subplots()
					x = data[:,0]
					y = data[:,1]
					ax.plot(x,y,'r')
					if args.xlim is not None:
						ax.set_xlim(args.xlim)
					if args.ylim is not None:
						ax.set_ylim(args.ylim)
					
					cout = stdt.eval(data,ax)
					
					if args.label is not None:
						for i,pp in enumerate(stdt.params):
							p = stdt.intp(i,cout[i])
							ind = np.floor(p[0]).astype(int)
							ax.text(x[ind],p[1],pp[0][0])

					fig.savefig(os.path.join(args.plotDir,os.path.splitext(os.path.basename(file))[0]+'.png'))
					plt.close(fig)	
				else:
					cout = stdt.eval(data,None)
			
				row = [os.path.basename(file)]+list(map(lambda x:'{:.2f}'.format(x),cout))
				print(fmtstr.format(*row))
				output.append(row)
		
		if args.fileout is not None:
			if args.fileout.endswith('.txt'):
				with open(args.fileout,'w',newline='') as f:
					f.write('\n'.join(['\t'.join(o) for o in output]))
			elif args.fileout.endswith('.csv'):
				with open(args.fileout,'w',newline='') as f:
					f.write('\n'.join([','.join(o) for o in output]))
			
	elif args.filein.endswith('.txt'):
		#single file mode
		lcdata = LCData(args.filein)
		data = lcdata.query(args.header)
		if args.noBLE:
			data[:,1] -= baselineMedian(data[:,1])
		else:
			data[:,1] -= baseline(data[:,1])
		
		x = data[:,0]
		y = data[:,1]
	
		fig,ax = plt.subplots()
		ax.plot(x,y,'r')
		if args.xlim is not None:
			ax.set_xlim(args.xlim)
		if args.ylim is not None:
			ax.set_ylim(args.ylim)
		cout=stdt.eval(data,ax)
		if args.label is not None:
			for i,pp in enumerate(stdt.params):
				p = stdt.intp(i,cout[i])
				ind = np.floor(p[0]).astype(int)
				ax.text(x[ind],p[1],pp[0][0])
		print()
		res = [['file']+[p[0][0] for p in stdt.params]]
		res.append([os.path.basename(args.filein)]+list(cout))
		print(tabulate(res,headers='firstrow'))
		plt.show()

