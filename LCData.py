import numpy as np
from tabulate import tabulate
import re
from Colors import Colors

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
	
	def getHeader(self):
		for c in self.content:
			if c.type == "LC Chromatogram":
				return c.name
		return None


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
		print(Colors.GREEN + self.name + Colors.RESET)
		#print(colors.RED + str(self.table) + colors.RESET)
		print(Colors.YELLOW + tabulate(self.table) + Colors.RESET)
		
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