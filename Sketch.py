# %config InlineBackend.figure_formats = ['svg']
import shapely.geometry as sg
from shapely.ops import unary_union
import matplotlib.pyplot as plt

class Draw:
	"""docstring for ClassName"""
	def __init__(self,balok,loc=(0,0)):
		self.balok = balok
		self.bm = balok.__dict__
		bm = balok.__dict__
		self.unit = bm["unit"]
		self.name = bm["name"]
		self.phi = bm["phi"]
		self.classification = bm["classification"]
		self.barloc = bm["barloc"]
		self.b = bm["bw"]
		self.h = bm["h"]
		self.ds = bm["ds"]
		self.b = bm["bw"]
		self.list_bar = bm["list_bar"]
		self.longBar = bm["longBar"]
		self.dRow = bm["dRow"]
		self.clearHSpace = bm["clearHSpace"]
		self.HSpacing = bm["HSpacing"]
		self.noDv = bm["noDv"]
		self.Dv = bm["Dv"]
		self.sv = bm["sv"]
		self.legs = bm["legs"]
		self.s = bm["s"]
		self.dt = bm["dt"]
		self.d = bm["d"]
		self.d1 = round(self.h-self.dt,2)
		self.a = bm["a"]
		self.c = bm["c"]
		self.es = bm["es"]
		self.et = bm["et"]
		self.e_cu = 0.003
		self.e_ty = 0.002
		self.be = 0
		self.hf = 0

		if hasattr(balok,"Cc"):
			self.Cc = bm["Cc"]
		if hasattr(balok,"Ccf"):
			self.Ccf = bm["Ccf"]
		if hasattr(balok,"Ccw"):
			self.Ccw = bm["Ccw"]

		self.T = bm["T"]

		if self.unit == "imperial":
			self.mm = 1/24.5
			self.unit_length = "in"
			self.unit_area = "in^2"
			self.unit_stress = "psi"
			self.unit_force = "kips"
			self.unit_moment = "kips-ft"
		else:
			self.mm = 1
			self.unit_length = "mm"
			self.unit_area = "mm^2"
			self.unit_stress = "MPa"
			self.unit_force = "kN"
			self.unit_moment = "kNm"

		if hasattr(balok,"be"):
			self.be = bm["be"]
			self.hf = bm["hf"]

		if hasattr(balok,"Dc"):
			self.Dc = bm["Dc"]
			self.nc = bm["nc"]
			self.dc = bm["dc"]
			self.Cs = bm["Cs"]
			self.es_c = bm["es_c"]
			self.noBarC = bm["noBarC"]

		self.noSlab = bm["noSlab"]
		self.nslab = bm["nslab"]
		self.Dslab = bm["Dslab"]

		self.add_d = 0*self.mm

	def boundarySmall(self,padX=300,padY=120):
		padX = padX*self.mm
		padY = padY*self.mm
		x0,y0 = self.loc
		x1,y1 = x0-(self.be/2)+(self.b/2)-padX, y0-padY/2
		x2,y2 = x0+self.b+(self.be/2)-(self.b/2)+padX, y0+self.h+padY
		x,y = sg.box(x1,y1,x2,y2).exterior.xy
		self.axs.fill(x, y, alpha=0., fc='white', ec='black')

	def boundaryLarge(self,padX=200,padY=220):
		padX = padX*self.mm
		padY = padY*self.mm
		padTop = 0.5*padY*self.mm
		x0,y0 = self.loc
		x_add = x0 + max(self.b,self.be/2) + 1.2*self.h + (1*self.h)
		x1,y1 = x0-(self.be/2)+(self.b/2)-padX, y0-padY
		x2,y2 = x_add+padX, y0+self.h+padTop
		x,y = sg.box(x1,y1,x2,y2).exterior.xy
		self.axs.fill(x, y, alpha=0.0, fc='white', ec='black')

	def beam(self,length):
	    x0,y0 = self.loc
	    x2 = x0+length
	    y2 = y0+self.h
	    x,y = sg.box(x0,y0,x2,y2).exterior.xy
	    self.axs.plot(x, y, color="black",linewidth=0.5)

	def dimV(self,dim,cordA,cordB,pad=50,rotate=0,alv="center",alh="right"):
		if self.hf != 0:
			if self.c < self.hf:
				alv = "bottom"

		start = (cordA[0]-pad*self.mm,cordA[1]-self.add_d)
		end = (cordB[0]-pad*self.mm,cordB[1]+self.add_d)
		textCoord = (start[0],start[1]+(end[1]-start[1])/2)
		plt.annotate(text='', xy=start, xytext=end, arrowprops=dict(arrowstyle='|-|',linewidth=0.5))
		plt.annotate(text='', xy=start, xytext=end, arrowprops=dict(arrowstyle='<|-|>',linewidth=0.5))
		plt.annotate(text=f'{dim} {self.unit_length}',xy=textCoord,rotation=rotate,va=alv,ha=alh,fontsize=self.fz)

	def dimH(self,dim,cordA,cordB,pad=50):
	    start = (cordA[0]-self.add_d,cordA[1]-pad*self.mm)
	    end = (cordB[0]+self.add_d,cordB[1]-pad*self.mm)
	    textCoord = (start[0]+(end[0]-start[0])/2,start[1])
	    plt.annotate(text='', xy=start, xytext=end, arrowprops=dict(arrowstyle='|-|',linewidth=0.5))
	    plt.annotate(text='', xy=start, xytext=end, arrowprops=dict(arrowstyle='<|-|>',linewidth=0.5))
	    plt.annotate(text=f'{dim} {self.unit_length}',xy=textCoord,rotation=0,va='bottom',ha="center",fontsize=self.fz)

	def rectangular(self):
	    # concrete section
	    x0,y0 = self.loc
	    sect01 = sg.box(x0,y0,x0+self.b,y0+self.h)
	    x,y = sect01.exterior.xy
	    self.axs.fill(x, y, alpha=1.0, fc='lightgrey', ec='black')

	def TBeam(self,dimension=True):
		x0,y0 = self.loc
		sect01 = sg.box(x0,y0,x0+self.b,y0+self.h)
		x1, y1 = x0-(self.be/2)+(self.b/2), self.h-self.hf
		x2, y2 = x0+self.b+(self.be/2)-(self.b/2), self.h
		sect02 = sg.box(x1,y1,x2,y2)
		merged = unary_union([sect01,sect02])
		x,y = merged.exterior.xy
		self.axs.fill(x, y, alpha=1.0, fc='lightgrey', ec='black')
		if self.unit == "imperial":
			pad = 1500*self.mm
		else:
			pad = 80*self.mm
		if dimension == True:
			self.dimV(f"{self.h}",(x1,y0),(x1,y2),pad)
			self.dimV(f"{self.hf}",(x1+self.be,y0+self.h-self.hf),(x1+self.be,y2),-pad)
			self.dimH(f"{self.b}",(x0,y0),(x0+self.b,y0),pad)
			self.dimH(f"{self.be}",(x1,y2),(x1+self.be,y2),-pad)

	def stirrups(self,notation=True):
		# stirrups
		x0,y0 = self.loc
		n_side = round(self.n/2)
		if self.be > 500*self.mm:
			x1 = x0+self.ds
			x_text = x0-self.h/6
			align = "right"
		else:
			x1 = x0+self.b-self.ds
			x_text = x0+self.b+self.h/6
			align = "left"
		y1 = y0+self.ds+self.h/2
		if notation == True:
			plt.annotate(text=f'{self.noDv}-{self.s}, {self.legs} legs', xy=(x1,y1),ha=align, xytext=(x_text,y1),color="black",fontsize=self.fz, arrowprops=dict(arrowstyle='->',color="black",linewidth=0.5))
			plt.annotate(text=f'Cover = {self.ds} {self.unit_length}', xy=(x1,y1-self.ds), ha=align, xytext=(x_text,y1-(2*self.ds)),color="black",fontsize=self.fz, arrowprops=dict(arrowstyle='->',color="black",linewidth=0.5))
		for leg in range(self.legs-1):
			if leg < n_side:
				sect02 = sg.box(x0+self.d1+(leg*self.sH),y0+self.d1,x0+self.b-self.d1,y0+self.h-self.d1)
			else:
				sect02 = sg.box(x0+self.d1,y0+self.d1,x0+self.d1+(leg*self.sH),y0+self.h-self.d1)
			x,y = sect02.buffer(self.Dv/2+self.D/2).exterior.xy
			self.axs.plot(x, y,color="black",linewidth=1.5)

	def dimension(self):
	    # dimension
	    pad = 80
	    x0,y0 = self.loc
	    x1,y1 = (x0,y0+self.h)
	    self.dimV(f"{self.h}",(x0,y0),(x1,y1),pad) # height

	    x2,y2 = (x0+self.b,y0)
	    self.dimH(f"{self.b}",(x0,y0),(x2,y2),pad) # width

	def barH(self,n,s,yrow,D,noBar,notation=True):
		x0, y0 = self.loc
		if self.barloc == "top":
			y0 = y0+self.h
			yrow = -yrow
		for i in range(0,n):
			shi = i*s
			bar = sg.Point(self.d1+shi, y0+yrow).buffer(D/2)
			# bar = sg.Point(x0+self.d1+shi, y0+y).buffer(D/2)
			x, y = bar.exterior.xy
			self.axs.fill(x, y, alpha=1, fc='grey', ec='black')
		pad = 10*self.mm
		x1 = x0+self.b-self.ds
		y1 = y0+yrow
		if notation == True:
			plt.annotate(text=f'{n}{noBar}', xy=(x1,y1), xytext=(x1+self.b/3,y1),va="bottom",ha="left", color="black",fontsize=self.fz, arrowprops=dict(arrowstyle='->',color="black",linewidth=0.5))

	def T_barH(self,n,s,yrow,D,noBar,notation=True):
		x0, y0 = self.loc
		if self.barloc == "top":
			y0 = y0+self.h
			yrow = -yrow
		for i in range(0,n):
			shi = i*s
			bar = sg.Point(self.d1+shi, y0+yrow).buffer(D/2)
			# bar = sg.Point(x0+self.d1+shi, y0+y).buffer(D/2)
			x, y = bar.exterior.xy
			self.axs.fill(x, y, alpha=1, fc='grey', ec='black')
		pad = 10*self.mm
		x1 = x0+self.b-self.ds
		y1 = y0+yrow
		if notation==True:
			plt.annotate(text=f'{n}{noBar}', xy=(x1,y1), xytext=(x1+self.b/3,y1-self.b/2.5), color="black",va="center",fontsize=self.fz, arrowprops=dict(arrowstyle='->',color="black",linewidth=0.5))

	def barHSlab(self,sh_slab,yrow,notation=True):
		x0,y0 = self.loc
		if self.barloc == "top":
			y0 = y0+self.h
			yrow = -yrow
		n_side = round(self.nslab/2)
		x1 = x0 + self.d1 - (n_side*sh_slab)
		x2 = x0 + self.b - self.d1 + sh_slab
		for i in range(0,n_side):
			shi = i*sh_slab
			bar = sg.Point(x1+shi, y0+yrow).buffer(self.Dslab/2)
			x, y = bar.exterior.xy
			self.axs.fill(x, y, alpha=1, fc='grey', ec='black')
		if notation==True:
			plt.annotate(text=f'{n_side}{self.noSlab}', xy=(x1,y0+yrow), xytext=(x1-(2*self.hf),y0+yrow-self.hf-self.d1), color="black",fontsize=self.fz, arrowprops=dict(arrowstyle='->',color="black",linewidth=0.5))
		for i in range(0,n_side):
			shi = i*sh_slab
			bar = sg.Point(x2+shi, y0+yrow).buffer(self.Dslab/2)
			x, y = bar.exterior.xy
			self.axs.fill(x, y, alpha=1, fc='grey', ec='black')
		if notation==True:
			plt.annotate(text=f'{n_side}{self.noSlab}', xy=(x2+((n_side-1)*sh_slab),y0+yrow), xytext=(x2+((n_side-1)*sh_slab)+(self.hf),y0+yrow-self.hf-self.d1), color="black",fontsize=self.fz, arrowprops=dict(arrowstyle='->',color="black",linewidth=0.5))

	def makeLine(self,start,end):
		plt.annotate(text='', xy=end, xytext=start, color="black")

	def diagrams(self):
		x0,y0 = self.loc
		x_label, ylabel = self.loc

		if self.barloc == "top":
			y0 = y0+self.h
			neg = -1
		else:
			neg = 1

		x_centerStrain = x0 + max(self.b,self.be/2) + 1.2*self.h

		coordA = (x_centerStrain, y0+(neg*(self.h-self.c)))
		coordB = (x_centerStrain, y0+(neg*self.h))

		self.dimV(f"c =\n{self.c}",coordA,coordB,60)

		if self.et > 0.005:
			strain_info = r"$\epsilon_{t} > 0.005$"
		elif self.et > 0.002:
			strain_info = r"$\epsilon_{ty}=0.002 < \epsilon_{t} < 0.005$"
		else:
			strain_info = r"$\epsilon_{t} < \epsilon_{ty}=0.002$"

		self.axs.text(x_centerStrain, ylabel-self.h/6, f"{strain_info}\n{self.classification}", rotation=0,ha="center", va='top',fontsize = self.fz)
		self.axs.text(x_centerStrain, ylabel-self.h/2.5, f"Strain Dist.", rotation=0,ha="center", va='top',fontsize = self.fz, fontweight="bold")

		# strain scaling
		ecu = 0.003
		scale = (self.b/2)/self.et
		x_ecu = x_centerStrain + (ecu*scale)
		x_et = x_centerStrain - (self.et*scale)

		poly = sg.Polygon([[x_centerStrain,y0+(neg*self.d1)], [x_et,y0+(neg*self.d1)], [x_ecu,y0+(neg*self.h)], [x_centerStrain,y0+(neg*self.h)],[x_centerStrain,y0]])
		x, y = poly.exterior.xy
		self.axs.plot(x, y, color='blue',linewidth=0.5)	# poly strain

		self.axs.plot((x_centerStrain,x_centerStrain),(y0,y0+(neg*self.h)),linewidth=1,color="black") # center strain
		self.axs.plot((x_centerStrain-self.b/2,x_centerStrain+self.b/2),(y0+(neg*(self.h-self.c)),y0+(neg*(self.h-self.c))),linewidth=1,color="black",linestyle="dashdot")
		et = r"$\epsilon_{t}$"
		plt.annotate(text=f'{et} = {round(self.et,4)}', xy=(x_centerStrain,y0), xytext=(x_centerStrain-self.ds,y0),ha="right", va="center", color="black",fontsize=self.fz)
		e_cu = r"$\epsilon_{cu}$"
		plt.annotate(text=f'{e_cu} = {ecu}', xy=(x_ecu,y0+(neg*(self.h+self.d1/4))), xytext=(x_ecu+(self.b/10),y0+(neg*self.h+self.d1/4)),ha="center",va="center", color="black",fontsize=self.fz)

		# force scaling
		x_centerForce = x_centerStrain + (1*self.h)
		self.axs.text(x_centerForce, ylabel-self.h/2.5, f"Int. Forces", rotation=0,ha="center", va='top',fontsize = self.fz, fontweight="bold")

		coordA = (x_centerForce, y0+(neg*(self.h-self.d)))
		coordB = (x_centerForce, y0+(neg*(self.h)))
		self.dimV(f"d =\n{self.d}",coordA,coordB,60,alv="center",alh="right")

		T = round(self.T/1000,2)

		scale = (self.b)/T
		x_T = x_centerForce + (T*scale)
		self.axs.plot((x_centerForce,x_centerForce),(y0,y0+(neg*self.h)),linewidth=1,color="black") # center strain
		plt.annotate(text=f'T = {T}', xy=(x_centerForce,y0+(neg*(self.h-self.d))), xytext=(x_T,y0+(neg*(self.h-self.d))),ha="left",va="center", color="black",fontsize=self.fz,arrowprops=dict(arrowstyle='-|>',color="red",linewidth=1))

		if hasattr(self,"Cc"):
			Cc = round(self.Cc/1000,1)
			x_Cc = x_centerForce + (Cc*scale)
			plt.annotate(text=f'Cc = {Cc}', xy=(x_centerForce,y0+(neg*(self.h-self.a/2))), xytext=(x_Cc,y0+(neg*(self.h-self.a/2))),ha="left",va="center", color="black",fontsize=self.fz,arrowprops=dict(arrowstyle='<|-',color="green",linewidth=1))

		if hasattr(self,"Cs"):
			Cs = round(self.Cs/1000,1)
			x_Cs = x_centerForce + (Cs*scale)
			plt.annotate(text=f'Cs = {Cs}', xy=(x_centerForce,y0+(neg*(self.h-self.dc))), xytext=(x_Cs,y0+(neg*(self.h-self.dc))),ha="left",va="center", color="black",fontsize=self.fz,arrowprops=dict(arrowstyle='<|-',color="green",linewidth=1))

		if hasattr(self,"Ccw"):
			Ccw = round(self.Ccw/1000,1)
			Ccf = round(self.Ccf/1000,1)
			x_Ccw = x_centerForce + (Ccw*scale)
			x_Ccf = x_centerForce + (Ccf*scale)
			plt.annotate(text=f'Ccw = {Ccw}', xy=(x_centerForce,y0+(neg*(self.h-self.a/2))), xytext=(x_Ccw,y0+(neg*(self.h-self.a/2))),ha="left",va="center", color="black",fontsize=self.fz,arrowprops=dict(arrowstyle='<|-',color="green",linewidth=1))
			plt.annotate(text=f'Ccf = {Ccf}', xy=(x_centerForce,y0+(neg*(self.h-self.hf/2))), xytext=(x_Ccf,y0+(neg*(self.h-self.hf/2))),ha="left",va="center", color="black",fontsize=self.fz,arrowprops=dict(arrowstyle='<|-',color="green",linewidth=1))

	def info(self):
		x0,y0=self.loc
		x_start = x0+max(self.b,self.be)/1.5 + 600*self.mm
		y_start = y0 + self.h - 10*self.mm

		# boudaries
		x,y = sg.box(x_start,y_start,x_start+500*self.mm,y0).exterior.xy
		self.axs.fill(x, y, alpha=0.0, fc='white', ec='black',linewidth=0.1)

		# Content
		self.axs.text(x_start,y_start, "Result:",va='center',fontsize = self.fz,fontweight='bold')

		y1, y2 =y_start-70*self.mm, y_start-70*self.mm
		x1, x2 =x_start, x_start+900*self.mm
		mg = 70*self.mm

		self.axs.text(x1,y1-(1*mg), f"Name = {self.bm['name']}",va='center',fontsize = self.fz)
		self.axs.text(x1,y1-(2*mg), f"Std = {self.bm['code']}",va='center',fontsize = self.fz)
		self.axs.text(x1,y1-(3*mg), f"Units = {self.bm['unit']}",va='center',fontsize = self.fz)
		self.axs.text(x1,y1-(4*mg), f"$f_c = {self.bm['fc']}\ {self.unit_stress}$",va='center',fontsize = self.fz)
		self.axs.text(x1,y1-(5*mg), f"$f_y = {self.bm['fy']}\ {self.unit_stress}$",va='center',fontsize = self.fz)
		self.axs.text(x1,y1-(6*mg), f"$E_s = {self.bm['Es']}\ {self.unit_stress}$",va='center',fontsize = self.fz)

		self.axs.text(x1,y2-(7*mg), f"$Classif. = {self.bm['classification']}$",va='center',fontsize = self.fz)
		self.axs.text(x1,y2-(8*mg), f"$Red.factor, \phi = {self.bm['phi']}$",va='center',fontsize = self.fz)
		self.axs.text(x1,y2-(9*mg), "Control:",va='center',fontsize = self.fz)
		if hasattr(self.balok,"Mu"):
			self.axs.text(x1+mg,y2-(10*mg), f"$M_u = {self.bm['Mu']}\ {self.unit_moment} < \phi M_n = {self.bm['phiMn']}\ {self.unit_moment}..OK$",va='center',fontsize = self.fz)
		else:
			self.axs.text(x1+mg,y2-(10*mg), f"$\phi M_n = {self.bm['phiMn']}\ {self.unit_moment}$",va='center',fontsize = self.fz)

		if hasattr(self.balok,"Vu"):
			self.axs.text(x1+mg,y2-(11*mg), f"$V_u = {self.bm['Vu']}\ {self.unit_force} < \phi V_n = {self.bm['phi_Vn']}\ {self.unit_force}..OK$",va='center',fontsize = self.fz)
		elif hasattr(self.balok,"phi_Vn"):
			self.axs.text(x1+mg,y2-(11*mg), f"$\phi V_n = {self.bm['phi_Vn']}\ {self.unit_force}$",va='center',fontsize = self.fz)

	def plot(self,loc=(0,0),size=(5,5),font=10,stirrups=False,diagram=False,supressinfo=False):
		self.fz = font
		self.loc = loc
		self.fig, self.axs = plt.subplots(figsize=size)
		self.axs.set_aspect('equal', 'datalim')
		if diagram == True:
			self.boundaryLarge()
		else:
			self.boundarySmall()
		# plot section title
		x0, y0 = self.loc
		self.axs.text(x0+self.b/2, y0-self.h/2.5, f"Cross Section\n{self.name}", rotation=0,ha="center", va='top',fontsize = self.fz, fontweight="bold")

		if supressinfo == False:
			if diagram == True:
				self.diagrams()
			else:
				self.info()

			if self.be == 0:
				self.rectangular()
				self.dimension()
				for i, row in enumerate(self.longBar):
					noBar, n = row
					D = self.list_bar[noBar][0]
					y = self.dRow[i]
					s = self.HSpacing[i]
					self.barH(n,s,y,D,noBar)
			else:
				self.TBeam()
				for i, row in enumerate(self.longBar):
					noBar, n = row
					D = self.list_bar[noBar][0]
					y = self.dRow[i]
					s = self.HSpacing[i]
					self.T_barH(n,s,y,D,noBar)

			# plot tension reinforcement bar
			noBar = self.longBar[0][0]
			self.n = self.longBar[0][1]
			self.D = self.list_bar[noBar][0]
			self.sH = self.HSpacing[0]

			if stirrups == True:
				self.stirrups()

			if hasattr(self,"Dc"):
				sc = (self.b-(2*(self.Dv+self.ds+self.Dc/2)))/ (self.nc-1)
				yrow = self.h-self.dc
				self.barH(self.nc,sc,yrow,self.Dc,self.noBarC)

			if self.nslab != 0:
				if self.barloc == "top":
					sh_slab = 150*self.mm
					yrow = self.h - self.dt
					self.barHSlab(sh_slab,yrow)
		else:
			if self.be == 0:
				self.rectangular()
				for i, row in enumerate(self.longBar):
					noBar, n = row
					D = self.list_bar[noBar][0]
					y = self.dRow[i]
					s = self.HSpacing[i]
					self.barH(n,s,y,D,noBar,notation=False)
			else:
				self.TBeam(dimension=False)
				for i, row in enumerate(self.longBar):
					noBar, n = row
					D = self.list_bar[noBar][0]
					y = self.dRow[i]
					s = self.HSpacing[i]
					self.T_barH(n,s,y,D,noBar,notation=False)

			# plot tension reinforcement bar
			noBar = self.longBar[0][0]
			self.n = self.longBar[0][1]
			self.D = self.list_bar[noBar][0]
			self.sH = self.HSpacing[0]

			if stirrups == True:
				self.stirrups(notation=False)

			if hasattr(self,"Dc"):
				sc = (self.b-(2*(self.Dv+self.ds+self.Dc/2)))/ (self.nc-1)
				yrow = self.h-self.dc
				self.barH(self.nc,sc,yrow,self.Dc,self.noBarC,notation=False)

			if self.nslab != 0:
				if self.barloc == "top":
					sh_slab = 150*self.mm
					yrow = self.h - self.dt
					self.barHSlab(sh_slab,yrow,notation=False)



		self.ax = plt.gca()
		self.ax.set_aspect('equal', adjustable='box')
		self.ax.axis("off")
		plt.show()
