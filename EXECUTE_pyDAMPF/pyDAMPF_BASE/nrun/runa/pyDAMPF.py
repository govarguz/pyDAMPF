#!/usr/bin/env python3
#=======================pyDAMPF==============================
#The code written below is part of the pyDAMPF project and has
#been developed as a joint project between the IJS and UMSA. This 
#software will enable students and researchers from the 
#AFM community to easily plan AFM experiments.
#======================CREDITS==============================
#
# PI Horacio V. Guzman Ph.D.
# B.Sc. Willy Menacho N.
#
#We thank you in advance for sending your feedback and/or 
#suggestions to:
#             horacio.v.g@gmail.com
#
#======================pyDAMPF LICENSE==============================

#Copyright (C) 2022  Horacio V. Guzman and Willy Menacho N.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# What we observe is not nature itself, but nature exposed to
# our method of questioning. 
#                                  W. Heisenberg
############################################################
#====================================================================
######for start timer#####
import time
inicio = time.time()

#========================Llajtamasis========================
# IMPORT MODULES
import os
import numpy as np
import matplotlib.pyplot as plt

# Define the path to the OS, for running python
#ROOT = os.path.realpath(os.path.dirname(__file__))
# first we run the simulation

import mypyDAMPF
#==========Define programs constants================
pi=3.141592654
kb=1.380658e-23
mu0=4*pi*1.e-7
nmax=10

# pi=mypyDAMPF.mypmmoda.pi
# kb=mypyDAMPF.mypmmoda.kb
# mu0=mypyDAMPF.mypmmoda.mu0
# nmax=mypyDAMPF.mypmmoda.nmax
#==========Initial problem conditions===============
x=0
t=0
rk4order=2
#==========Initialize output vectors================
y=np.zeros(nmax)
yout=np.zeros(nmax)
#|||||||||||||||||||||||||||||||||||||||||||||||||||||
#=================INPUTS========================

#------------------AUTO GATHERING--------------------------
case = open("case.in", "r")
cinp = np.genfromtxt(case)
sid = cinp[6]
# same the same
omega0a = cinp[5]*1000.
omegaa = cinp[5]*1000.
q = cinp[1]
kc = cinp[4]
a0 = cinp[3]*1e-9
emuestra = cinp[2]*1000000.
epunta = 150.0*1e9
rad = cinp[0]*1e-9
a00 = 0.165*1e-9
ham = cinp[9]*1e-20
eta = cinp[8]
epsilon1 = cinp[16]	#60.0
epsilon2 = cinp[7]*1e-3
delta = 0.0
ljmin = cinp[22]*1e-9	#0.0
length = cinp[23]*1e-9	#0.0
mtip = 0.0
msample = 0.0
nper = cinp[10]	#1000.0
npp = cinp[11]		#1024.0
naux = int(nper*npp)
nperfin = cinp[12]	#128.0
dzc = cinp[14]*1e-9	#0.01e-9
if a0 < 2.1*1e-9:
    zcmax = a0*(1.0+1.5)
else:
    zcmax = a0*(1.0+0.1)
    
#zcmin = rad*(-1.0)
zcmin = cinp[13]*1e-9	#0
if emuestra > 0 and epunta > 0:        # NEW
    ebarra=(4.0/3.0)*emuestra
else:
    ebarra=0.
    
fixzc = cinp[15]*1e-9	#0.35*a0

#dlvo
eps= cinp[17]
eps0= cinp[18]*1e-12	#8.85e-12
sigmat= cinp[19]*1e-3	#0.032
sigmas= cinp[20]*1e-3	#0.05
debye= cinp[21]*1e-9	#0.48e-9
#=========================== Saving Sim Inputs ===============================
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
tx00=[ '"res freq 1 (kHz)"', '"drive freq 2 (kHz)"', '"k1 (N/m)"', '"Q1 (adim)"', '"A01 (V)"', '"Tip Radius (nm)"', '"a0 (nm)"', '"Reduced YM (Pa)"', '"Hamaker Const(J)"', '"Viscosity (Pa s)"','"LR Disp. (adim)"', '"Cap Disp.(adim)"', '"LJdepth (nN)"', '"LJlength (nm)"', '"Tip magn (A m2)"', '"Surf magn(A m2)"','"N periods osc"','"N point period"', '"N per to the ss"', '"zc max (nm)"', '"zc min (nm)"', '"delta zc (nm)"', '"zc fixed td(nm)"' ,'"rel perm  (C^2/N m^2)"','"vac perm (C^2/N m^2)"','"t-s density (mC/m^2)"','"s char den  (mC/m^2)"','"Debye (nm)"']
tx00=np.array(tx00)
tx01=[ omega0a, omegaa, kc, q, a0, rad, a00, ebarra, ham, eta, epsilon1, epsilon2, ljmin, length, mtip, msample, nper, npp, nperfin, zcmax, zcmin, dzc, fixzc ,eps,eps0,sigmat,sigmas,debye]
tx01=np.array(tx01)
f1 = open("inputs.ini","w")
for i in np.ndindex(tx00.shape):
	filerow='%s      %5.5e\n' % (tx00[i],tx01[i])
	f1.write(filerow)
f1.close()
#    elif: ADD FILE HAS BEEN LOADED IN ADVANCE
#=========================Call fortran functions ===============================
#print (naux)
mypyDAMPF.mypmmoda.mainbim(naux)
#print ("Runnig f90 module PM")
#=========================Files handling w Python ===============================
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

#palo2="dforcekara"
#os.mkdir(palo1)
sid=int(sid*1)
#print ("Sikita michaway Le Chardenal")
piat2=os.getcwd()
naind="inputs"+str(sid)+".ini"
nazcd="zcdom"+str(sid)+".dfo"
natd="tdom"+str(sid)+".dfo"
os.rename("inputs.ini",naind)
os.rename("zcdom.dfo",nazcd)
os.rename("tdom.dfo",natd)
#os.mkdir(palo1)
#========================Obtaining outputs zc domain ==========================
#========Reading files for point-mass outputs===========
#||||||||||||||||||||||||||||||||||||||||||
tdom = np.genfromtxt(natd)
tt = tdom[:,0]
zt = tdom[:,1]
vt = tdom[:,2]
forcet = tdom[:,4]
maxforcet = tdom[:,5]
vdwt = tdom[:,6]
hertzt = tdom[:,7]
viscot = tdom[:,8]
capt = tdom [:,9]
zcdom = np.genfromtxt(nazcd)
zc = zcdom[:,0]
amp1d = zcdom[:,1]
fase1d = zcdom[:,2]
dmind = zcdom[:,3]
dmaxd = zcdom[:,4] # new
fmaxd = zcdom[:,5]
defld = zcdom[:,9]
ets1d = zcdom[:,6]
vts1d = zcdom[:,7]
pts1d = zcdom[:,11] # new
fmedd = zcdom[:,8] # new
tcd = zcdom[:,10]# new
#==================================================
# End of getting arrays from simulated DATA
#======================PLOT========================
#========Building the big Switch===========
#||||||||||||||||||||||||||||||||||||||||||
'''
#def amp1():
#========Amp 1 vs. zc=================
plt.figure(1)
plt.grid(True)
plt.title('$A_1$($z_c$)')
plt.ylabel('$A_1 \; [nm]$',weight='bold',size='x-large')
plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
plt.plot(zc,amp1d,'r.--')
plt.xlim(0.5,10)
plt.savefig('amp1Vzc'+str(sid)+'.png')
plt.close()
#def pha1():
    #========Phase 1 vs. zc=================
plt.figure(2)
plt.grid(True)
plt.title('$\phi_1$($z_c$)')
plt.ylabel('$\phi_1 \; [deg]$',weight='bold',size='x-large')
plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
plt.plot(zc,fase1d,'b.--')
plt.xlim(0.5,10)
plt.savefig('phi1Vzc'+str(sid)+'.png')
plt.close()
'''
'''
#def den1():
    #========Diss Energy 1st vs. zc=================
plt.figure(3)
plt.grid(True)
plt.title('First disp. energy vs. Average distance')
plt.ylabel('$E_{ts1} * 10^{-20} \; [J]$',weight='bold',size='x-large')
plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
plt.plot(zc,ets1d,'g.--')
plt.xlim(0.5,10)
plt.savefig('ets1Vzc'+str(sid)+'.png')
plt.close()
#def dpo1():
    #========Diss Power 1st vs. zc=================
plt.figure(4)
plt.grid(True)
plt.title('First disp. power vs. Average distance')
plt.ylabel('$P_{ts1} * 10^{-20} \; [W]$',weight='bold',size='x-large')
plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
plt.plot(zc,pts1d,'r.--')
plt.xlim(0.5,10)
plt.savefig('pts1Vzc'+str(sid)+'.png')
plt.close()
#def vir1():
    #========Virial 1st vs. zc=================
plt.figure(5)
plt.grid(True)
plt.title('Virial first vs. Average distance')
plt.ylabel('$V_{ts1} * 10^{-20} \; [J]$',weight='bold',size='x-large')
plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
plt.plot(zc,vts1d,'b.--')
plt.xlim(0.5,10)
plt.savefig('vts1Vzc'+str(sid)+'.png')
plt.close()
#def def1():
    #========Deflection vs. zc=================
plt.figure(6)
plt.grid(True)
plt.title('Deflection vs. Average distance')
plt.ylabel('$Defl \; [nm]$',weight='bold',size='x-large')
plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
plt.plot(zc,defld,'r.--')
plt.xlim(0.5,10)
plt.savefig('deflVzc'+str(sid)+'.png')
plt.close()
#def dmin1():
    #========Dmin vs. zc=================
plt.figure(7)
plt.grid(True)
plt.title('Minimum Distance vs. Average distance')
plt.ylabel('$dmin \; [nm]$',weight='bold',size='x-large')
plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
plt.plot(zc,dmind,'g.--')
plt.xlim(0.5,10)
plt.savefig('dminVzc'+str(sid)+'.png')
plt.close()
#def maf():
    #========Fmax vs. zc=================
plt.figure(8)
plt.grid(True)
plt.title('Maximum Force vs. Average distance')
plt.ylabel('$F_{max} \; [nN]$',weight='bold',size='x-large')
plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
plt.plot(zc,fmaxd,'b.--')
plt.xlim(0.5,10)
plt.savefig('fmaxVzc'+str(sid)+'.png')
plt.close()
#def tco():
    #========Contact time vs. zc=================
plt.figure(9)
plt.grid(True)
plt.title('Contact Time (tc/T) vs. Average distance (nm)')
plt.ylabel('$Contact Time \; [t/t_{pixel\;imaging}]$',weight='bold',size='x-large')
plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
plt.plot(zc,tcd,'r.--')
plt.xlim(0.5,10)
plt.savefig('ctimeVzc'+str(sid)+'.png')
plt.close()
'''
'''
def amp2():
        #========Amp 2 vs. zc=================
    plt.figure(10)
    plt.grid(True)            
    plt.title('$A_2$($z_c$)')
    plt.ylabel('$A_2 \; [nm]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    plt.plot(zc,amp2d,'y.--')
    plt.savefig('amp2Vzc'+str(sid)+'.png')
def pha2():
    #========Phase 2 vs. zc=================
    plt.figure(11)
    plt.grid(True)    
    plt.title('$\phi_2$($z_c$)')
    plt.ylabel('$\phi_2 \; [deg]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    plt.plot(zc,fase2d,'g.--')
    plt.savefig('phi2Vzc'+str(sid)+'.png')
def den2():
    #========Diss Energy 2nd vs. zc=================
    plt.figure(12)
    plt.grid(True)
    plt.title('Second disp. energy vs. Average distance')
    plt.ylabel('$E_{ts2} * 10^{-20} \; [J]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    plt.plot(zc,ets2d,'r.--')
    plt.savefig('ets2Vzc'+str(sid)+'.png')
def dpo2():
    #========Diss Power 1st vs. zc=================
    plt.figure(13)
    plt.grid(True)
    plt.title('Second disp. power vs. Average distance')
    plt.ylabel('$P_{ts2} * 10^{-20} \; [W]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    plt.plot(zc,pts2d,'r.--')
    plt.savefig('pts2Vzc'+str(sid)+'.png')
def vir2():
    #========Virial 2nd vs. zc=================
    plt.figure(14)
    plt.grid(True)        
    plt.title('Virial second vs. Average distance ')
    plt.ylabel('$V_{ts2} * 10^{-20} \; [J]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    plt.plot(zc,vts2d,'g.--')
    plt.savefig('vts2Vzc'+str(sid)+'.png')
'''
    #----------------------------------------------------------------
    #===================Time Domain==================================
    #----------------------------------------------------------------
'''
#def ztt():
    #========Instantaneous position vs. time=================
plt.figure(15)
plt.grid(True)
plt.title("$z_T$($t/T$)")
plt.ylabel('$z_T \; [nm]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,zt,'b.--')
plt.savefig('totalposTime'+str(sid)+'.png')
plt.close()
#def vtt():
    #========Instantaneous velocity vs. time=================
plt.figure(16)
plt.grid(True)
plt.title("$v_T$($t/T$)")
plt.ylabel('$v_T \; [nm/\mu s]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,vt,'y.--')
plt.savefig('totalveloTime'+str(sid)+'.png')
plt.close()
#def ftt():
    #========Instantaneous force vs. time=================
plt.figure(17)
plt.grid(True)
plt.title('$f_T$($t/T$)')
plt.ylabel('$f_T \; [nN]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,forcet,'r.--')
plt.savefig('totalforceTime'+str(sid)+'.png')
plt.close()
#def fttz():
    #==========Instantaneous force vs. time (zoom)========
plt.figure(18)
plt.grid(True)
plt.title('$f_T$($t/T$)')
plt.ylabel('$f_T \; [nN]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,forcet,'r.--')
plt.xlim(nper-3,nper)
plt.savefig('totalforceTimezoom'+str(sid)+'.png')
plt.close()
#def ftdi():
    #=========max force vs. distance (ida)===============
plt.figure(19)
aux=np.where(forcet[0:1024]==forcet[0:1024].max())
aux=int(aux[0][0])
plt.grid(True)
plt.title('$f max_G$($d$)')
plt.ylabel('$f max \; [nN]$',weight='bold',size='x-large')
plt.xlabel('$[d]$',weight='bold',size='x-large')
plt.plot(maxforcet[0:aux],forcet[0:aux],'b.--')
plt.savefig('totalmaxforceTimeida'+str(sid)+'.png')
plt.close()
#def ftdv():
    #=========max force vs. distance (vuelta)===============
plt.figure(20)
aux=np.where(forcet[0:1024]==forcet[0:1024].max())
aux=int(aux[0][0])
plt.grid(True)
plt.title('$f max_R$($d$)')
plt.ylabel('$f max \; [nN]$',weight='bold',size='x-large')
plt.xlabel('$[d]$',weight='bold',size='x-large')
plt.plot(maxforcet[aux:1024],forcet[aux:1024],'g.--')
plt.savefig('totalmaxforceTimevuelta'+str(sid)+'.png')
plt.close()
#def ftvdw():
    #=========vdw force vs. time============================
plt.figure(21)
plt.grid(True)
plt.title('$f vdw$($d$)')
plt.ylabel('$f vdw \; [nN]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,vdwt,'b.--')
plt.xlim(nper-3, nper)
plt.savefig('vdwforce'+str(sid)+'.png')
plt.close()
#def fthertz():
    #=========hertz force vs. time=========================
plt.figure(22)
plt.grid(True)
plt.title('$f hertz$($d$)')
plt.ylabel('$f hertz \; [nN]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,hertzt,'b.--')
plt.xlim(nper-3, nper)
plt.savefig('hertzforce'+str(sid)+'.png')
plt.close()
#def ftvisco():
    #=========viscoelastic force vs. time==================
plt.figure(23)
plt.grid(True)
plt.title('$f viscoelastica$($d$)')
plt.ylabel('$f visco \; [nN]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,viscot,'b.--')
plt.xlim(nper-3, nper)
plt.savefig('viscoforce'+str(sid)+'.png')
plt.close()
#def ftcap():
    #=========capilar force vs. time============================
plt.figure(24)
plt.grid(True)
plt.title('$f capilar$($d$)')
plt.ylabel('$f cap \; [nN]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,capt,'b.--')
plt.xlim(nper-3, nper)
plt.savefig('capforce'+str(sid)+'.png')
plt.close()
#def ftdt():
    #=========max force vs. distance (ida-vuelta)===============
plt.figure(25)
plt.grid(True)
plt.title('$f max TOTAL$($d$)')
plt.ylabel('$f max \; [nN]$',weight='bold',size='x-large')
plt.xlabel('$[d]$',weight='bold',size='x-large')
plt.plot(maxforcet[0:1024],forcet[0:1024],'c.--')
plt.savefig('totalmaxforceTimeidavuelta'+str(sid)+'.png')
plt.close()

#def zttzoom():
    #========Instantaneous position vs. time(zoom)=================
plt.figure(26)
plt.grid(True)
plt.title("$z_T$($t/T$)")
plt.ylabel('$z_T \; [nm]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,zt,'b.--')
plt.xlim(nper-3, nper)
plt.savefig('totalposTimezoom'+str(sid)+'.png')
plt.close()
#def vttzoom():
    #========Instantaneous velocity vs. time(zoom)=================
plt.figure(27)
plt.grid(True)
plt.title("$v_T$($t/T$)")
plt.ylabel('$v_T \; [nm/\mu s]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,vt,'y.--')
plt.xlim(nper-3, nper)
plt.savefig('totalveloTimezoom'+str(sid)+'.png')
plt.close()
#def fttzoom():
    #========Instantaneous force vs. time(zoom)=================
plt.figure(28)
plt.grid(True)
plt.title('$f_T$($t/T$)')
plt.ylabel('$f_T \; [nN]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,forcet,'r.--')
plt.xlim(nper-3, nper)
plt.savefig('totalforceTimezoom'+str(sid)+'.png')
plt.close()

#def p1a1():
    #========Phase 1 vs. Amp 21=================
plt.figure(29)
plt.grid(True)
plt.title("$\phi_1$($A_1$)")
plt.ylabel('$\phi_1 \; [deg]$',weight='bold',size='x-large')
plt.xlabel('$A_1\; [nm]$',weight='bold',size='x-large')
plt.plot(amp1d,fase1d,'bo')
plt.savefig('phi1Vamp1'+str(sid)+'.png')
plt.close()


#def zt1():
    #========Instantaneous position 1 vs. time=================
plt.figure(30)
plt.grid(True)
plt.title("$z_1$($t/T$)")
plt.ylabel('$z_1 \; [nm]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,zt,'b.--')
plt.savefig('FirstposTime'+str(sid)+'.png')
plt.close()
#def zt2():
    #========Instantaneous position 2 vs. time=================
plt.figure(31)
plt.grid(True)
plt.title("$z_2$($t/T$)")
plt.ylabel('$z_2 \; [nm]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,zt,'g.--')
plt.savefig('SecondposTime'+str(sid)+'.png')
plt.close()
#def zt3():
    #========Instantaneous position 3 vs. time=================
plt.figure(32)
plt.grid(True)
plt.title("$z_3$($t/T$)")
plt.ylabel('$z_3 \; [nm]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,zt,'y.--')
plt.savefig('ThirdposTime'+str(sid)+'.png')
plt.close()
#def zt4():
    #========Instantaneous position 4 vs. time=================
plt.figure(33)
plt.grid(True)
plt.title("$z_4$($t/T$)")
plt.ylabel('$z_4 \; [nm]$',weight='bold',size='x-large')
plt.xlabel('$[t/T]$',weight='bold',size='x-large')
plt.plot(tt,zt,'g.--')
plt.savefig('FourthposTime'+str(sid)+'.png')
plt.close()

def p2a2():
    #========Phase 2 vs. Amp 2=================
    plt.figure(34)
    plt.grid(True)
    plt.title("$\phi_2$($A_2$)")
    plt.ylabel('$\phi_2 \; [deg]$',weight='bold',size='x-large')
    plt.xlabel('$A_2\; [nm]$',weight='bold',size='x-large')
    plt.plot(amp2d,fase2d,'ro')
    plt.savefig('phi2Vamp2'+str(sid)+'.png')
    plt.close()
def p2a1():
    #========Phase 2 vs. Amp 1=================
    plt.figure(35)
    plt.grid(True)
    plt.title("$\phi_2$($A_1$)")
    plt.ylabel('$\phi_2 \; [deg]$',weight='bold',size='x-large')
    plt.xlabel('$A_1\; [nm]$',weight='bold',size='x-large')
    plt.plot(amp1d,fase2d,'ro')
    plt.savefig('phi2Vamp1'+str(sid)+'.png')
    plt.close()
def p1a2():
    #========Phase 1 vs. Amp 2=================
    plt.figure(36)
    plt.grid(True)
    plt.title("$\phi_1$($A_2$)")
    plt.ylabel('$\phi_1 \; [deg]$',weight='bold',size='x-large')
    plt.xlabel('$A_2\; [nm]$',weight='bold',size='x-large')
    plt.plot(amp2d,fase1d,'go')
    plt.savefig('phi1Vamp2'+str(sid)+'.png')
    plt.close()
   
#Create a Tuple for the command glade_ui.get_widget
#curves = [amp1, pha1, den1, dpo1, vir1, def1, dmin1, maf, tco, amp2, pha2, den2, dpo2, vir2, ztt, vtt, ftt, zt1, zt2, zt3, zt4, p1a1, p2a2, p2a1, p1a2]
#opt1=[1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
#||||||||||||||||||||||||||||||||||||||||||
#=========Big switch ready to use==============
print ("This is Pujllay")
for ii in range(len(opt1)):
    if opt1[ii]==1:
        curves[ii]()    #calling ccurves.py and saving plots
    else:
        print ("You did not visualize", ii+1)    #-----Message
'''  

#####for end timer#####
fin = time.time()
#print('EXECUTION TIME CASE '+str(sid)+': ', (fin-inicio)/60 , '[min]') 
time = (fin-inicio)/60
file = open('time'+str(sid)+'.txt','w')
file.write('time (seg)=    %s'%time)
file.close()
#print('FINISH CASE'+str(sid))    
