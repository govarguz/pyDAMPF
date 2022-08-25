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
'''
######for start timer#####
import time
inicio = time.time()
'''
# IMPORT MODULES
import os
import shutil
# IMPORT NEWLY VISUALIZATION MODULES
import numpy as np

#======================Gathering=============================

base = np.genfromtxt('tempall.txt')
iden = base[:,0]
f0base = base[:,4]
rbase = base[:,5]
kbase = base[:,6]
nperbase = base[:,9]
nppbase = base[:,10]
nperfinbase = base[:,11]
zcminbase = base[:,12]
dzcbase = base[:,13]
fixzcbase = base[:,14]
a0base = base[:,15]
ebase = base[:,16]
qbase = base[:,17]
hambase = base[:,18]
etabase = base[:,19]
epsilon2base = base[:,20]
epsilon1base = base[:,21]
epsbase = base[:,22]
eps0base = base[:,23]
sigmatbase = base[:,24]
sigmasbase = base[:,25]
debyebase = base[:,26]
ljminbase = base[:,27]
lenghtbase = base[:,28]

tcases = np.size(iden)

#======================Parallelization=============================

factor=1  
paraleliz=int(tcases/factor)        

#------------------Case inputs gen--------------------------
for i in range(0,paraleliz):
    # placeholders for progress checking
    idname='case'+str(i+1)+'.in'
    np.savetxt(idname, [0])

for i in range(0,paraleliz):
    txu=np.array([rbase[i],qbase[i],ebase[i],a0base[i],kbase[i],f0base[i],iden[i],epsilon2base[i],etabase[i],hambase[i],nperbase[i],nppbase[i],nperfinbase[i],zcminbase[i],dzcbase[i],fixzcbase[i],epsilon1base[i],epsbase[i],eps0base[i],sigmatbase[i],sigmasbase[i],debyebase[i],ljminbase[i],lenghtbase[i]])
    txu=txu.T
    idname='case'+str(i+1)+'.in'
    np.savetxt(idname,txu)
    # after saving
    piat=os.getcwd()
    piat1=piat
    dirname='runa'
    baname='case.in'
    shutil.copy(idname,piat+'/'+dirname)
    os.chdir(piat+'/'+dirname)
    os.rename(idname,baname)
    #os.remove(idname)
    # text has been saved now is time to run the code
    print('will run pydampf')
    exec(open("pyDAMPF.py").read())
    print('did run pydampf')
    #exec(compile(open("dforcenv1.py").read(), "dforcenv1.py", 'exec'))
    #%run dforcenv1.py
    os.chdir(piat1)
    
#####for end timer#####
'''
fin = time.time()
print('EXECUTION TIME THREAD '+str(factor)+': ', (fin-inicio)/60 , '[min]') 
time1 = (fin-inicio)/60
file = open('timeTHREAD'+str(factor)+'.txt','w')
file.write('time=%s'%time1)
file.close()
'''
#print('FINISH THREAD  '+str(factor))  


