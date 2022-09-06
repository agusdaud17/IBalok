"""
This code for design and analyzing reinforced conrete beam.
By: Agus Daud
"""


import Reinforcement
import SNI2847
import ACI318
from sympy import*
from math import*

class Rect:
	"""docstring for """
	def __init__(self,name="B1",*kwargs):
		self.material()
		self.name = name

	def material(self,fc=30, fy=400, Es=200000, code="SNI 2847:2019"):
		self.fc = fc
		self.fy = fy
		self.Es = Es
		self.code = code

		if self.code == "SNI 2847:2019":
			self.std = SNI2847
			self.unit = "metric"
			self.list_bar = Reinforcement.metric_bar
			self.mm = 1
			self.to_kNm = 1/1000000
		else:
			self.std = ACI318
			self.unit = "imperial"
			self.list_bar = Reinforcement.imperial_bar
			self.mm = 1/24.5
			self.to_kNm = 1/12/1000

	def dimension(self,bw=300,h=None,ds=40,decimal=2,barloc="bottom"):
		self.bw = bw
		self.h = h
		self.ds = ds
		self.decimal = decimal
		self.barloc = barloc

	def rebar(self,sv=25,Dv="D10",legs=2,nslab=0,Dslab="D13",s=150,**kwargs):
		if "d" in kwargs:
			self.d = kwargs['d']
			self.dt = kwargs['d']
			self.effective_height = "by user"
		else:
			self.effective_height = "by program"
		if "As" in kwargs:
			self.As = kwargs['As']

		if "rebar" in kwargs:
			self.longBar = kwargs['rebar']
			self.row = len(self.longBar)

		if self.unit == "imperial":
			Dv = "#3"

		self.noDv = Dv
		self.Dv = self.list_bar[Dv][0]
		self.legs = legs
		self.sv = round(sv*self.mm,1) # clear vertical spacing
		self.s = round(s*self.mm)

		self.nslab = nslab
		if self.nslab == 0:
			if self.unit == "imperial":
				Dslab = "#4"
		self.noSlab = Dslab
		self.Dslab, self.Abar_slab = self.list_bar[Dslab]
		self.As_slab = round(self.nslab*self.Abar_slab,self.decimal)

	def effectiveDepth(self,b):
		self.AsRow = []
		self.NBar = []
		self.clearHSpace = []
		self.HSpacing = []
		self.dRow = []
		self.nRow = []
		self.DRow = []

		for row in self.longBar:
			numBar = row[1]
			noBar = row[0]
			D = self.list_bar[noBar][0]
			self.rearrange(numBar,noBar,D)

		self.longBar = self.nRow

		MRow = 0

		for i, row in enumerate(self.longBar):
			noBar = row[0]
			numBar = row[1]
			D, As = self.list_bar[noBar]
			As_tot = As*numBar
			hs_clear = round((self.bw - (2*(self.ds+self.Dv)) - (numBar*D))/(numBar-1),self.decimal)
			hs = round((self.bw - (2*(self.ds+self.Dv+D/2)))/(numBar-1),self.decimal)

			if i == 0:
				d = self.ds + self.Dv + D/2
			else:
				d = self.dRow[i-1] + max(self.sv,D) + D

			self.AsRow.append(As_tot)
			self.dRow.append(d)
			self.NBar.append(numBar)
			self.clearHSpace.append(hs_clear)
			self.HSpacing.append(hs)
			MRow += As_tot*d

		self.As = round(sum(self.AsRow),self.decimal)

		g = MRow/self.As

		if self.nslab != 0:
			self.As = round(sum(self.AsRow) + self.As_slab,self.decimal)
		else:
			self.As = round(sum(self.AsRow),self.decimal)

		dt = self.h - self.dRow[0] 		# distance from compression face to the bar closes to tension face.
		self.d = round(self.h - g,self.decimal)	# Over all effective depth.
		self.Ag = round(self.bw*self.d,self.decimal)
		self.rho = round(self.As/self.Ag,6)
		return dt

	def compressBlockDepth(self,b):
		a = round((self.As*self.fy) / (0.85*self.fc*b),self.decimal)
		return a

	@property
	def neutral(self):
		c = round(self.a/self.beta1,self.decimal)
		return c

	@property
	def tensileStrain(self):
		et = round(((self.dt-self.c)/self.c * self.std.e_cu),6)
		return et

	@property
	def tensileStrainCentroid(self):
		es = round(((self.d-self.c)/self.c * self.std.e_cu),6)
		return es

	def run(self):
		# if hasattr(self,'d'):
		# 	pass
		# else:
		# 	self.dt = self.effectiveDepth

		b = self.bw
		if hasattr(self,"dt"):
			pass
		else:
			self.dt = self.effectiveDepth(b=b)

		self.a = self.compressBlockDepth(b=b)

		self.beta1 = round(self.std.beta_1(self.fc),self.decimal)
		self.c = self.neutral
		self.es = self.tensileStrainCentroid

		if abs(self.es) <= self.std.e_ty:
			As, Es, d, b, c, ecu, fc, beta1 = symbols('A_s E_s d b c e_cu f_c beta1')
			As = self.As
			Es = self.Es
			d = self.d
			b = self.bw
			ecu = self.std.e_cu
			fc = self.fc
			beta1 = self.beta1

			T = As*Es*ecu* ((d-c) / c) # Page 138/155 pdf reference book
			C_c = 0.85*fc*b*beta1*c
			result = solveset(Eq(T-C_c,0),c)

			# recalculated
			self.c = round(max(result.__dict__['_args_set']),self.decimal)
			self.a = self.compressBlockDepth(b)
			self.Cc = round(0.85*self.fc*b*self.beta1*self.c,self.decimal)
			self.T = round(self.As*self.Es*self.std.e_cu*((self.d-self.c)/self.c),self.decimal)
			self.Mn = round(self.T*(self.d-self.a/2)*self.to_kNm,self.decimal)
			if self.effective_height == "by user":
				pass
			else:
				self.dt = self.effectiveDepth(b)
		else:
			self.Cc = round((0.85*self.fc*b*self.beta1*self.c),self.decimal)
			self.T = round((self.As*self.fy),self.decimal)
			self.Mn = round(self.T*(self.d-self.a/2)*self.to_kNm,self.decimal)

		self.et = self.tensileStrain
		self.phi, self.classification = self.std.strengthReductionFactor(self.et)
		self.phiMn = round(self.phi*self.Mn,self.decimal)
		# self.factoredMomentStrength()

	def checkingCondition(self):
		if self.clearHSpace[-1] > self.sv+(0*self.mm):
			# increase number of N bar
			self.longBar[-1][1] += 1
		else:
			# Increase number of row
			db = self.longBar[0][0]
			self.longBar.append([db,2])

	def controlReinfRatio(self):
		self.rho_min, self.et_max = self.std.standardBeam(self.fc,self.fy)
		if self.et < self.et_max:
			self.condition = "et < et_max, overreinforce!"
		elif self.rho > self.rho_min:
			self.condition = "rho > rho_min ..OK"
		else:
			self.condition = "rho_min > rho, revise beam reinforcement!"

	def rearrange(self,numBar,noBar,D):
		# check spacing
		lapis = []
		n = 2

		while True:
			n_sisa = numBar - sum(lapis)
			s_tul = (self.bw-(2*(self.ds+self.Dv))-(2*2*self.Dv)-(n*D)) / (n-1)
			if n_sisa == 0:
				break
			elif n_sisa ==1:
				self.nRow.append([noBar,n_sisa+1])
				break
			elif n_sisa < n:
				lapis.append(n_sisa)
				self.nRow.append([noBar,n_sisa])
				break
			if s_tul < max(self.sv,D):
				lapis.append(n)
				self.nRow.append([noBar,n])
				n = 1
			n += 1

	def design(self,Mu,D,Dv="D10"):
		# self.rebar()
		self.noDv = Dv
		self.Dv = self.list_bar[Dv][0]

		# Assume single layer of reinforcement page 214 pdf
		self.Mu = Mu
		noBar = D
		j = 0.95
		D, As_tul = self.list_bar[D]
		self.sv = max(25*self.mm,D)
		d = self.h - self.ds -self.Dv - D/2
		class_assumed = "Tension-controlled"
		phi = 0.9 # Assume tension controlled
		As = (self.Mu/self.to_kNm) / (phi*self.fy*j*d)
		numBar = max(ceil(As/As_tul),2)

		# check spacing
		self.nRow = []
		self.rearrange(numBar,noBar,D)

		self.longBar = self.nRow

		while True:
			self.rebar(rebar=self.longBar)
			self.run()

			if self.phiMn >= self.Mu:
				break
			else:
				if self.classification == "Transition":
					self.h += round(50*self.mm,0)
				else:
					self.checkingCondition()
				self.run()

	property
	def factoredMomentStrength(self):
		phi_Mn = self.phi*self.Mn
		self.phiMn = round(phi_Mn,self.decimal)
		return phi_Mn

	def shear(self,Vu,legs=2):
		self.legs = legs
		self.Vu = Vu
		self.phi_v, self.Vc, self.Vmax1, self.Vmax2, self.Avmin_s, self.Av_req, self.stirrup_req, self.section_req = self.std.shearStd(self.Vu,self.fc,self.bw,self.d,self.fy)

		if self.section_req != "OK":
			raise Exception("Sorry, section not large enough! increase beam depth (h)")

		self.Atul_v = round(0.25*3.14*self.Dv**2,self.decimal)
		self.s = round(250*self.mm)

		while True:
			self.Av_prov = round((self.legs*self.Atul_v)/self.s,self.decimal)
			self.Vs = round((self.Av_prov*self.fy*self.d)/1000,self.decimal)
			self.s_max = self.std.maxStirrupSpace(self.Vs,self.Vmax1,self.d)
			self.Vn = round(self.Vc+self.Vs,self.decimal)
			self.phi_Vn = round(self.phi_v*self.Vn,self.decimal)

			if self.Vu < self.phi_Vn:
				break
			else:
				if self.s > 100*self.mm:
					self.s -= round(25*self.mm)
				else:
					self.legs += 1
					self.s = round(250*self.mm)

class Doubly(Rect):
	"""docstring for Doubly"""
	def __init__(self,name="B1"):
		super().__init__()
		self.name = name

	def rebar(self,sv=25,Dv="D10",legs=2,nc=0,Dc=0,dsc=0,nslab=0,Dslab="D13",s=150,**kwargs):
		if "d" in kwargs:
			self.d = kwargs['d']
			self.dt = kwargs['d']
		if "As" in kwargs:
			self.As = kwargs['As']

		if "rebar" in kwargs:
			self.longBar = kwargs['rebar']
			self.row = len(self.longBar)

		self.Dc = self.list_bar[Dc][0]
		self.noBarC = Dc
		self.nc = nc
		self.Abar_c = self.list_bar[Dc][1]
		self.Asc = round(self.nc*self.Abar_c,self.decimal)

		if self.unit == "imperial":
			self.noDv = "#3"
		else:
			self.noDv = Dv
		self.Dv = self.list_bar[self.noDv][0]
		self.legs = legs
		self.sv = round(sv*self.mm,1) # clear vertical spacing
		self.dc = dsc+self.Dv+self.Dc/2
		self.s = round(s*self.mm)

		self.nslab = nslab
		if self.nslab ==0:
			if self.unit == "imperial":
				Dslab = "#4"

		self.noSlab = Dslab
		self.Dslab, self.Abar_slab = self.list_bar[Dslab]
		self.As_slab = round(self.nslab*self.Abar_slab,self.decimal)

	def run(self):
		b = self.bw
		self.dt = round(self.effectiveDepth(b=b),self.decimal)
		self.beta1 = round(self.std.beta_1(self.fc),self.decimal)

		# start solving quadrad equation
		# assuming tension steel yield and; compression steel not yield
		As, Es, d, b, c, ecu, fc, beta1, es_c = symbols('A_s E_s d b c e_cu f_c beta1 e_sc')
		As = self.As
		Es = self.Es
		d = self.d
		b = self.bw
		ecu = self.std.e_cu
		fc = self.fc
		beta1 = self.beta1

		es_c = ((c-self.dc)/c)*ecu
		fsc = Es*es_c
		T = As*self.fy
		C_c = 0.85*fc*b*beta1*c
		C_s = self.Asc*(fsc-0.85*fc)
		result = solveset(Eq(T-C_c-C_s,0),c)

		self.c = round(max(result.__dict__['_args_set']),self.decimal)
		self.es_c = round(((self.c-self.dc)/self.c)*self.std.e_cu,6)
		self.fsc = round(self.Es*self.es_c,self.decimal)
		self.es = self.tensileStrainCentroid

		# recalculated
		self.a = self.beta1*self.c
		self.Cc = round(0.85*self.fc*b*self.beta1*self.c,self.decimal)
		self.Cs = round(self.Asc*(self.fsc-(0.85*self.fc)),self.decimal)
		self.T = round(self.As*self.fy,self.decimal)

		if self.es < self.std.e_ty+self.std.e_cu:
			# start solving quadrad equation
			# tension steel not yiels
			As, Es, d, b, c, ecu, fc, beta1, es_c = symbols('A_s E_s d b c e_cu f_c beta1 e_sc')
			As = self.As
			Es = self.Es
			d = self.d
			b = self.bw
			ecu = self.std.e_cu
			fc = self.fc
			beta1 = self.beta1

			es_c = ((c-self.dc)/c)*ecu
			fsc = Es*es_c
			T = As*Es*ecu* ((d-c) / c)
			C_c = 0.85*fc*b*beta1*c
			C_s = self.Asc*(fsc-0.85*fc)
			result = solveset(Eq(T-C_c-C_s,0),c)

			self.c = round(max(result.__dict__['_args_set']),self.decimal)
			self.es_c = round(((self.c-self.dc)/self.c)*self.std.e_cu,6)
			self.fsc = round(self.Es*self.es_c,self.decimal)
			self.es = self.tensileStrainCentroid

			# recalculated
			self.a = self.beta1*self.c
			self.Cc = round(0.85*self.fc*b*self.beta1*self.c,self.decimal)
			self.Cs = round(self.Asc*(self.fsc-(0.85*self.fc)),self.decimal)
			self.T = round(self.As*self.Es*self.es,self.decimal)

			if self.fsc > self.fy:
				# tension steel not yiels; and compression steel yield
				As, Es, d, b, c, ecu, fc, beta1, es_c = symbols('A_s E_s d b c e_cu f_c beta1 e_sc')
				As = self.As
				Es = self.Es
				d = self.d
				b = self.bw
				ecu = self.std.e_cu
				fc = self.fc
				beta1 = self.beta1

				es_c = ((c-self.dc)/c)*ecu
				fsc = self.fy
				T = As*Es*ecu* ((d-c) / c)
				C_c = 0.85*fc*b*beta1*c
				C_s = self.Asc*(fsc-0.85*fc)
				result = solveset(Eq(T-C_c-C_s,0),c)

				self.c = round(max(result.__dict__['_args_set']),self.decimal)
				self.es_c = round(((self.c-self.dc)/self.c)*self.std.e_cu,6)
				self.fsc = round(self.fy,self.decimal)
				self.es = self.tensileStrainCentroid

				# recalculated
				self.a = self.beta1*self.c
				self.Cc = round(0.85*self.fc*b*self.beta1*self.c,self.decimal)
				self.Cs = round(self.Asc*(self.fsc-(0.85*self.fc)),self.decimal)
				self.T = round(self.As*self.Es*self.es,self.decimal)

		# self.Mn = round((self.Cc*(self.d-self.a/2)) + (self.Cs*(self.d-self.dc)) *self.to_kNm,self.decimal)
		self.dt = self.effectiveDepth(b)
		self.et = self.tensileStrain
		self.es = self.tensileStrainCentroid
		self.phi, self.classification = self.std.strengthReductionFactor(self.et)
		self.Mn = round((self.Cc*(self.d - self.a/2) + self.Cs*(self.d - self.dc))*self.to_kNm,self.decimal)
		self.phiMn = round(self.phi*self.Mn,self.decimal)

	def design(self,Mu,D,Dc,Dv="D10",percent=0.5):
		self.Dc, self.Abar_c = self.list_bar[Dc]
		self.noDv = Dv
		self.Dv = self.list_bar[Dv][0]
		# self.rebar()

		# Assume single layer of reinforcement page 214 pdf
		self.Mu = Mu
		noBar = D
		j = 0.95
		D, As_tul = self.list_bar[D]
		self.sv = max(25*self.mm,D)
		d = self.h - self.ds - D/2
		class_assumed = "Tension-controlled"
		phi = 0.9 # Assume tension controlled
		As = (self.Mu/self.to_kNm) / (phi*self.fy*j*d)
		numBar = ceil(As/As_tul)

		# check spacing
		self.nRow = []
		self.rearrange(numBar,noBar,D)

		n_c = ceil((percent*As)/self.Abar_c)
		self.nc = max(n_c,2)

		self.longBar = self.nRow

		while True:
			self.rebar(rebar=self.longBar,nc=self.nc,Dc=Dc)
			self.run()
			self.controlReinfRatio()

			if self.phiMn > self.Mu:
				break
			else:
				if self.classification == "Transition":
					self.h += round(50*self.mm,0)
				else:
					self.checkingCondition()
				self.run()

class Flanged(Rect):
	"""docstring for Flanged"""
	def __init__(self,name="B1"):
		super().__init__()
		self.name = name

	def dimension(self,bw=300,h=None,hf=0,be=0,ds=40,decimal=2,barloc="bottom"):
		self.bw = bw
		self.h = h
		self.hf = hf
		self.be = be
		self.ds = ds
		self.decimal = decimal
		self.legs = 2
		self.barloc = barloc

	def run(self):
		if self.barloc == "top":
			b = self.bw
			self.analysis = Rect.run(self)
			self.analysis
		else:
			b = self.be
			self.dt = self.effectiveDepth(b=b)
			self.a = self.compressBlockDepth(b=b)
			self.beta1 = round(self.std.beta_1(self.fc),self.decimal)
			self.c = self.neutral

			if self.a < self.hf:
				self.et = self.tensileStrain
				self.es = self.tensileStrainCentroid
				self.phi, self.classification = self.std.strengthReductionFactor(self.et)
				self.Cc = round(0.85*self.fc*b*self.a,self.decimal)
				self.T = round((self.As+self.As_slab)*self.fy,self.decimal)
				self.Mn = round(self.T*(self.d-self.a/2)*self.to_kNm,self.decimal)
				self.phiMn = round(self.phi*self.Mn,self.decimal)
			else:
				self.Ccf = 0.85*self.fc*(self.be-self.bw)*self.hf
				self.T = (self.As+self.As_slab)*self.fy
				self.a = round((self.T-self.Ccf)/(0.85*self.fc*self.bw),self.decimal)
				self.Ccw = 0.85*self.fc*self.bw*self.a
				self.c = self.c = self.neutral
				self.et = self.tensileStrain
				self.es = self.tensileStrainCentroid
				self.phi, self.classification = self.std.strengthReductionFactor(self.et)
				self.Mn = round((self.Ccf*(self.d-self.hf/2) + self.Ccw*(self.d-self.a/2))*self.to_kNm,self.decimal)
				self.phiMn = round(self.phi*self.Mn,self.decimal)

	def design(self,Mu,D,Dv="D10"):
		if self.unit == 'imperial':
			Dv = "#3"

		self.rebar(Dv=Dv)
		self.noDv = Dv
		self.Dv = self.list_bar[Dv][0]

		if self.barloc == "top":
			b = self.bw
		else:
			b = self.be

		# Assume single layer of reinforcement page 214 pdf
		self.Mu = Mu
		noBar = D
		j = 0.95
		D, As_tul = self.list_bar[D]
		self.sv = max(25*self.mm,D)
		d = self.h - self.ds - D/2
		class_assumed = "Tension-controlled"
		phi = 0.95 # Assume tension controlled
		As = (self.Mu/self.to_kNm) / (phi*self.fy*j*d)
		numBar = ceil(As/As_tul)

		# check spacing
		self.nRow = []
		self.rearrange(numBar,noBar,D)

		self.longBar = self.nRow

		while True:
			self.rebar(rebar=self.longBar,Dv=Dv)

			if self.barloc == "top":
				self.analysis = Rect.run(self)
				self.analysis
			else:
				self.run()

			self.controlReinfRatio()

			if self.phiMn > self.Mu:
				break
			else:
				if self.classification == "Transition":
					self.h += round(50*self.mm,0)
				else:
					self.checkingCondition()

				if self.barloc == "top":
					self.analysis
				else:
					self.run()

class DoublyFlanged(Doubly):
	"""docstring for Flanged"""
	def __init__(self,name="B1"):
		super().__init__()
		self.name = name

	def dimension(self,bw=300,h=None,hf=0,be=0,ds=40,decimal=2,barloc="bottom"):
		self.bw = bw
		self.h = h
		self.hf = hf
		self.be = be
		self.ds = ds
		self.decimal = decimal
		self.legs = 2
		self.barloc = barloc

	def run(self):
		if self.barloc == "top":
			b = self.bw
			self.dt = round(self.effectiveDepth(b=b),self.decimal)
			self.beta1 = round(self.std.beta_1(self.fc),self.decimal)

			# start solving quadrad equation
			As, Es, d, b, c, ecu, fc, beta1, es_c = symbols('A_s E_s d b c e_cu f_c beta1 e_sc')
			As = self.As+self.As_slab
			Es = self.Es
			d = self.dt
			b = self.bw
			ecu = self.std.e_cu
			fc = self.fc
			beta1 = self.beta1

			es_c = ((c-self.dc)/c)*ecu
			fsc = Es*es_c
			T = As*self.fy
			C_c = 0.85*fc*b*beta1*c
			C_s = self.Asc*(fsc-0.85*fc)
			result = solveset(Eq(T-C_c-C_s,0),c)

			# recalculated
			self.c = round(max(result.__dict__['_args_set']),self.decimal)
			self.es_c = round(((self.c-self.dc)/self.c)*self.std.e_cu,6)
			self.fsc = round(self.Es*self.es_c,self.decimal)
			if self.fsc >= self.fy:
				self.fsc = self.fy
			self.a = self.compressBlockDepth(b)
			self.Cc = round(0.85*self.fc*b*self.beta1*self.c,self.decimal)
			self.Cs = round(self.Asc*(self.fsc-(0.85*self.fc)),self.decimal)
			self.T = round(self.As*self.fy,self.decimal)
			self.Mn = round((self.Cc*(self.d-self.a/2)) + (self.Cs*(self.d-self.dc)) *self.to_kNm,self.decimal)
			self.dt = self.effectiveDepth(b)
			self.et = self.tensileStrain
			self.es = self.tensileStrainCentroid
			self.phi, self.classification = self.std.strengthReductionFactor(self.et)
			self.Mn = round((self.Cc*(self.d - self.a/2) + self.Cs*(self.d - self.dc))*self.to_kNm,self.decimal)
			self.phiMn = round(self.phi*self.Mn,self.decimal)
		else:
			b = self.be
			self.dt = self.effectiveDepth(b=b)
			self.a = self.compressBlockDepth(b=b)
			self.beta1 = round(self.std.beta_1(self.fc),self.decimal)
			self.c = self.neutral

			if self.a < self.hf:
				# ingnore compression reinforcement page 156 reference book
				self.es_c = 0
				self.fsc = 0
				self.Cc = round(0.85*self.fc*b*self.beta1*self.c,self.decimal)
				self.Cs = 0
				self.T = round(self.As*self.fy,self.decimal)
				self.et = self.tensileStrain
				self.es = self.tensileStrainCentroid
				self.phi, self.classification = self.std.strengthReductionFactor(self.et)
				self.Mn = round(self.T*(self.d - self.a/2)*self.to_kNm,self.decimal)
				self.phiMn = round(self.phi*self.Mn,self.decimal)
			else:
				# start solving quadrad equation
				As, Es, d, b, c, ecu, fc, beta1, es_c = symbols('A_s E_s d b c e_cu f_c beta1 e_sc')
				As = self.As
				Es = self.Es
				d = self.dt
				b = self.be
				ecu = self.std.e_cu
				fc = self.fc
				beta1 = self.beta1

				es_c = ((c-self.dc)/c)*ecu
				fsc = Es*es_c
				T = As*self.fy
				C_s = self.Asc*(fsc-0.85*fc)
				Ccw = 0.85*self.fc*self.bw*self.beta1*c
				Ccf = 0.85*self.fc*(self.be-self.bw)*self.hf
				result = solveset(Eq(T-Ccw-Ccw-C_s,0),c)

				# recalculated
				self.c = round(max(result.__dict__['_args_set']),self.decimal)
				self.es_c = round(((self.c-self.dc)/self.c)*self.std.e_cu,6)
				self.fsc = round(self.Es*self.es_c,self.decimal)
				self.a = self.compressBlockDepth(b)
				self.Ccw = round(0.85*self.fc*self.bw*self.beta1*self.c,self.decimal)
				self.Ccf = round(0.85*self.fc*(self.be-self.bw)*self.hf,self.decimal)
				self.Cs = round(self.Asc*(self.fsc-(0.85*self.fc)),self.decimal)
				self.T = round(self.As*self.fy,self.decimal)
				self.dt = self.effectiveDepth(b)
				self.et = self.tensileStrain
				self.es = self.tensileStrainCentroid
				self.phi, self.classification = self.std.strengthReductionFactor(self.et)
				self.Mn = round((self.Ccf*(self.d-self.hf/2) + self.Ccw*(self.d-self.a/2) + self.Cs*(self.d-self.dc))*self.to_kNm,self.decimal)
				self.phiMn = round(self.phi*self.Mn,self.decimal)

	def design(self,Mu,D,Dc,Dv="D10",percent=0.5):
		# self.rebar()
		if self.barloc == "top":
			b = self.bw
		else:
			b = self.be
		self.Dc, self.Abar_c = self.list_bar[Dc]
		self.noDv = Dv
		self.Dv = self.list_bar[Dv][0]
		# self.rebar()

		# Assume single layer of reinforcement page 214 pdf
		self.Mu = Mu
		noBar = D
		j = 0.95
		D, As_tul = self.list_bar[D]
		self.sv = max(25*self.mm,D)
		d = self.h - self.ds - D/2
		class_assumed = "Tension-controlled"
		phi = 0.9 # Assume tension controlled
		As = (self.Mu/self.to_kNm) / (phi*self.fy*j*d)
		numBar = ceil(As/As_tul)

		# check spacing
		self.nRow = []
		self.rearrange(numBar,noBar,D)

		n_c = ceil((percent*As)/self.Abar_c)
		self.nc = max(n_c,2)

		self.longBar = self.nRow

		while True:
			self.rebar(rebar=self.longBar,nc=self.nc,Dc=Dc)
			self.run()
			self.controlReinfRatio()

			if self.phiMn > self.Mu:
				break
			else:
				if self.classification == "Transition":
					self.h += round(50*self.mm,0)
				else:
					self.checkingCondition()
				self.run()
