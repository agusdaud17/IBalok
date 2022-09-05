"""
Nilai berikut merupakan rekapitulasi peraturan beton berdasarkan SNI 1728:2019
Oleh: Agus Daud
Tgl: 13 Juli 2022
"""
from math import*

# Elastisitas beton
def Ec(fc):
	return 4700*sqrt(fc)

# Regangan tarik baja tulangan ulir, 21.2.2.1
# untuk fy 420 MPa, sama dengan fy/Es
e_ty = 0.002
e_cu = 0.003

# Faktor yang berhubungan tinggi blok tegangan ekivalen, 22.2.2.4.3
def beta_1(fc):
	if fc <= 28:
		beta1 = 0.85
	elif fc < 55:
		beta1 = 0.85 - (0.05*((fc-28)/7))
	else:
		beta1 = 0.65
	return beta1


# Faktor reduksi kekuatan (phi) Tabel 21.2.2
def strengthReductionFactor(et):
	if abs(et)<= e_ty:
		classification = "Compression-controlled"
		phi = 0.65
	elif abs(et) < e_ty+e_cu:
		classification = "Transition"
		phi = round(0.65 + (0.25*(abs(et)-e_ty)/(0.005-e_ty)),6)
	else:
		classification = "Tension-controlled"
		phi = 0.9
	return phi, classification

phi_shear = 0.75
phi_torsion = 0.75
phi_bearing = 0.65

# Standar tulangan longitudinal untuk balok
def standardBeam(fc,fy):
	rho_min = round(max(0.25*sqrt(fc)/fy, 1.4/fy),6)
	et_max = 0.004
	return rho_min, et_max

def shearStd(Vu,fc,b,d,fy):
	Vc = round(0.17*sqrt(fc)*b*d /1000,2)
	Vmax1 = round(0.33*sqrt(fc)*b*d /1000,2)
	Vmax2 = round(0.66*sqrt(fc)*b*d /1000,2)
	Avmin1 = round(0.062*sqrt(fc)*(b/fy),2)
	Avmin2 = round((0.35*b)/fy,2)
	Avmin_s = round(max(Avmin1,Avmin2),2)

	# Is the cross section large enough?
	if Vu < phi_shear*(Vc+Vmax2):
		section_req = "OK"
	else:
		section_req = "not large enough"

	# Are stirrups reequired by ACI/SNI code...?
	if Vu/phi_shear >= Vc:
		stirrup_req = "stirrups required"
		Av_req = (Vu - phi_shear*Vc)*1000 / (phi_shear*fy*d)
		Av_req = round(max(Av_req,Avmin_s),2)
	elif Vu/phi_shear > Vc/2:
		stirrup_req = "minimum stirrups provided"
		Av_req = round(Avmin_s,2)
	else:
		stirrup_req = "stirrups not required"
		Av_req = round(Avmin_s,2)

	return phi_shear, Vc, Vmax1, Vmax2, Avmin_s, Av_req, stirrup_req, section_req

def maxStirrupSpace(Vs,Vmax1,d):
	# max spacing
	if Vs <= Vmax1:
		s_max = min(d/2,600)
	else:
		s_max = min(d/4,300)
	return round(s_max,0)
