import argparse
import os
import re
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

from ParseList import parseList
from Standard import Standard
import Baseline
from LCData import LCData
from Colors import Colors

os.system('')

HEADER_A1 = "LC Chromatogram(Detector A-Ch1)"
HEADER_A2 = "LC Chromatogram(Detector A-Ch2)"
HEADER_B1 = "LC Chromatogram(Detector B-Ch1)"

def tryint(s):
	try:
		return int(s)
	except:
		return s

def alnum_key(s):
	return [tryint(c) for c in re.split('([0-9]+)',s)]

def loadDataBC(file,args):	
	data = loadData(file,args)
	if args.noBaselineCorrection:
		data[:,1] -= Baseline.baselineMedian(data[:,1])
	else:
		data[:,1] -= Baseline.baseline(data[:,1])
	return data

def loadData(file,args):
	lcdata = LCData(file)
	#if args.header is None:
	#	args.header = lcdata.getHeader()
	#	if args.header:
	#		print('header : ' + Colors.YELLOW + args.header + Colors.RESET)
	#	else:
	#		print(file + ' has no LC Chromatogram Data')
	#		return None

	data = lcdata.query(args.header)
	
	if args.polarity:
		data[:,1] *= -1
	
	n = len(data)
	tmax = data[-1,0]

	if not np.isclose(tmax/(n-1),1/120):
		nn = tmax * 120 + 1
		if not args.dataSizeChangedMessage:
			args.dataSizeChangedMessage = True
			print('data size changed from ' +  str(n) + ' to ' + str(nn))
		n = len(data)
		x = np.interp(np.linspace(0,n,nn),np.arange(n),data[:,0])
		y = np.interp(np.linspace(0,n,nn),np.arange(n),data[:,1])
		data = np.transpose([x,y])

	return data

parser = argparse.ArgumentParser(description='LC processor')
parser.add_argument('filein',nargs='?',help='Input file or directory')
parser.add_argument('-o','--fileout',help='Output file')
parser.add_argument('-s','--standardFile',help="Standard file")
parser.add_argument('-p','--paramFile',help="Parameter file")
parser.add_argument('--plotDir',help="Plot Directory")
parser.add_argument('--plotParamsDir',help='Plot parameters Directory')
parser.add_argument('--checkParams',action='store_true',help='check parameters (Standard file is required)')
parser.add_argument('--checkStandards',action='store_true',help='check standards data (Standard file is required)')
parser.add_argument('--header',default=HEADER_B1,help='header name to query in a datafile')
parser.add_argument('--noBaselineCorrection',action='store_true',help='no baseline correction')
parser.add_argument('--xlim',help='x region of plot')
parser.add_argument('--ylim',help='y region of plot')
parser.add_argument('--label',help='show chemical species label',action='store_true')
parser.add_argument('--peakProminence',type=float,default=1.0,help='set prominence parameter for peak detection. default : 1.0')
parser.add_argument('--peakWidth',type=float,default=0.1,help='set width parameter for peak detection. default : 0.1')
parser.add_argument('--peakInclude',help='set regions to include for peak detection. signature : [[tmin1,tmax2],[tmin1,tmax2],...]')
parser.add_argument('--peakExclude',help='set regions to exclude for peak detection. signature : [[tmin1,tmax2],[tmin1,tmax2],...]')
parser.add_argument('--polarity',action='store_true',help='multiply -1 to the data')
parser.add_argument('--shiftTolerance',type=float,help='set x shift tolerance')
parser.add_argument('--shift',default=1,type=float,help='set x shift multiplier')
parser.add_argument('--noTotalFit',action='store_true',help='no total fit')
parser.add_argument('--gaussianFit',action='store_true',help='fit with multiple gaussins')

args = parser.parse_args()
args.dataSizeChangedMessage = False

if args.header == 'A1':
	args.header = HEADER_A1
if args.header == 'A2':
	args.header = HEADER_A2

if args.xlim is not None:
	args.xlim = parseList(args.xlim)
if args.ylim is not None:
	args.ylim = parseList(args.ylim)

if args.peakInclude is not None:
	args.peakInclude = parseList(args.peakInclude)
	if args.peakInclude and not isinstance(args.peakInclude[0],list):
		args.peakInclude = [args.peakInclude]

if args.peakExclude is not None: 
	args.peakExclude = parseList(args.peakExclude)
	if args.peakExclude and not isinstance(args.peakExclude[0],list):
		args.peakExclude = [args.peakExclude]

stdt = Standard(args)

if args.standardFile is not None:
	stdt.load(args.standardFile,lambda file:loadData(file,args))
	stdt.checkParams()
	if args.paramFile is not None:
		stdt.saveParams()
	if args.checkStandards:
		stdt.checkStandards()
		exit()
	if not args.filein and args.plotDir is not None:
		stdt.plotParams()
		exit()
elif args.paramFile is not None:
	stdt.loadParams()
	if args.checkParams:
		stdt.checkParams()
		exit()
	if not args.filein and args.plotDir is not None:
		stdt.plotParams()
		exit()
elif not args.filein:
	parser.print_help()
	exit()
elif os.path.isdir(args.filein):
	files = os.listdir(args.filein)
	files.sort(key=alnum_key)
	fig,ax = plt.subplots()
	
	if args.xlim is not None:
		ax.set_xlim(args.xlim)
	if args.ylim is not None:
		ax.set_ylim(args.ylim)
	
	for file in files:
		if file.endswith('.txt'):
			data = loadDataBC(os.path.join(args.filein,file),args)
			x = data[:,0]
			y = data[:,1]
			ax.plot(x,y,label=os.path.basename(file),alpha=0.5)
	
	ax.legend()
	if args.plotDir is not None:
		os.makedirs(args.plotDir,exist_ok=True)
		fig.savefig(os.path.join(args.plotDir,os.path.splitext(os.path.basename(args.filein))[0]+'.png'))
		plt.close(fig)
	else:
		plt.show()
	exit()
else:
	data = loadDataBC(args.filein,args)
	x = data[:,0]
	y = data[:,1]

	fig,ax = plt.subplots()

	if args.xlim is not None:
		ax.set_xlim(args.xlim)
	if args.ylim is not None:
		ax.set_ylim(args.ylim)

	ax.plot(x,y)
	
	if args.plotDir is not None:
		os.makedirs(args.plotDir,exist_ok=True)
		fig.savefig(os.path.join(args.plotDir,os.path.splitext(os.path.basename(args.filein))[0]+'.png'))
		plt.close(fig)
	else:	
		plt.show()
	exit()	

if args.filein is not None:
	if os.path.isdir(args.filein):
		#multiple file mode	
		files = os.listdir(args.filein)
		files.sort(key=alnum_key)
		row = ['File Name']+list(stdt.chems.keys())
		fmtstr = "{:<20}"
		for chem in stdt.chems:
			fmtstr += " {:<15}"
		print(fmtstr.format(*row))
		output = []
		output.append(row)
		if args.plotDir is not None:
			os.makedirs(args.plotDir,exist_ok=True)
		for file in files:
			if file.endswith('.txt'):
				data = loadDataBC(os.path.join(args.filein,file),args)
			
				if args.plotDir is not None:
					fig,ax = plt.subplots()
					x = data[:,0]
					y = data[:,1]

					if args.xlim is not None:
						ax.set_xlim(args.xlim)
					if args.ylim is not None:
						ax.set_ylim(args.ylim)
					
					ax.plot(x,y)
					
					ax.set_title(os.path.basename(file))
					cout = stdt.eval(data,ax)
					
					if args.label is not None:
						for i,chem in enumerate(stdt.chems.values()):
							p = chem.interpolate(cout[i])
							ax.text(x[np.floor(chem.pos).astype(int)],p[1],chem.name)

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
		data = loadDataBC(args.filein,args)
		x = data[:,0]
		y = data[:,1]
	
		fig,ax = plt.subplots()
		ax.plot(x,y)
		if args.xlim is not None:
			ax.set_xlim(args.xlim)
		if args.ylim is not None:
			ax.set_ylim(args.ylim)
		cout=stdt.eval(data,ax)
		if args.label is not None:
			for i,chem in enumerate(stdt.chems.values()):
				p = chem.interpolate(cout[i])
				ax.text(x[np.floor(chem.pos).astype(int)],p[1],chem.name)
		print()
		res = [['file']+list(stdt.chems.keys())]
		res.append([os.path.basename(args.filein)]+list(cout))
		print(tabulate(res,headers='firstrow'))

		if args.plotDir is not None:
			fig.savefig(os.path.join(args.plotDir,os.path.splitext(os.path.basename(args.filein))[0]+'.png'))
			plt.close(fig)
		else:	
			plt.show()
