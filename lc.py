import argparse
import os
#import matplotlib.patches as patches

from Standard import *
from ParseList import parseList
from LCData import LCData

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

def loadData(file,args):
	data = LCData(file).query(args.header)
	if args.noBaselineCorrection:
		data[:,1] -= baselineMedian(data[:,1])
	else:
		data[:,1] -= baseline(data[:,1])
	return data

parser = argparse.ArgumentParser(description='LC processor')
parser.add_argument('filein',nargs='?',help='Input file or directory')
parser.add_argument('-o','--fileout',help='Output file')
parser.add_argument('-s','--standardFile',help="Standard file")
parser.add_argument('-p','--paramFile',help="Parameter file")
parser.add_argument('--plotDir',help="Plot Directory")
parser.add_argument('--plotParamsDir',help='Plot parameters Directory')
parser.add_argument('--checkParams',action='store_const',const=10,help='check parameters (Standard file is required)')
parser.add_argument('--checkStandards',action='store_true',help='check standards data (Standard file is required)')
parser.add_argument('--header',default=HEADER_B1,help='header name to query in a datafile')
parser.add_argument('--noBaselineCorrection',default=False,action='store_true',help='no baseline correction')
parser.add_argument('--xlim',help='x region of plot')
parser.add_argument('--ylim',help='y region of plot')
parser.add_argument('--label',help='show chemical species label',action='store_true')
parser.add_argument('--peakProminence',type=float,default=1.0,help='set prominence parameter for peak detection')
parser.add_argument('--peakWidth',type=int,default=15,help='set width parameter for peak detection')
parser.add_argument('--peakInclude',help='set regions to include for peak detection. signature : [[tmin1,tmax2],[tmin1,tmax2],...]')
parser.add_argument('--peakExclude',help='set regions to exclude for peak detection. signature : [[tmin1,tmax2],[tmin1,tmax2],...]')

args = parser.parse_args()

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
	stdt.load(args.standardFile,lambda file:LCData(file).query(args.header))
	print(Colors.YELLOW + "Standard Data Table" + Colors.RESET)
	print(tabulate(stdt.table))
	if args.paramFile is not None:
		stdt.saveParams()
	if args.checkParams:
		stdt.checkParams()
		exit()
	if args.checkStandards:
		stdt.checkStandards()
		exit()
	if args.plotParamsDir is not None:
		stdt.plotParams()
		exit()
elif args.paramFile is not None:
	stdt.loadParams()
	if args.checkParams:
		stdt.checkParams()
	if args.plotParamsDir is not None:
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
			data = loadData(os.path.join(args.filein,file),args)
			x = data[:,0]
			y = data[:,1]
			ax.plot(x,y,label=os.path.basename(file))
	
	ax.legend()
	if args.plotDir is not None:
		os.makedirs(args.plotDir,exist_ok=True)
		fig.savefig(os.path.join(args.plotDir,os.path.splitext(os.path.basename(args.filein))[0]+'.png'))
		plt.close(fig)
	else:	
		plt.show()
	exit()
else:
	data = loadData(args.filein,args)
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
				data = loadData(os.path.join(args.filein,file),args)
			
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
						for i,pp in enumerate(stdt.params):
							p = stdt.interpolate(i,cout[i])
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
		data = loadData(args.filein,args)
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
			for i,pp in enumerate(stdt.params):
				p = stdt.interpolate(i,cout[i])
				ind = np.floor(p[0]).astype(int)
				ax.text(x[ind],p[1],pp[0][0])
		print()
		res = [['file']+[p[0][0] for p in stdt.params]]
		res.append([os.path.basename(args.filein)]+list(cout))
		print(tabulate(res,headers='firstrow'))

		if args.plotDir is not None:
			fig.savefig(os.path.join(args.plotDir,os.path.splitext(os.path.basename(args.filein))[0]+'.png'))
			plt.close(fig)
		else:	
			plt.show()
