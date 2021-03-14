#!/usr/bin/env python3.7
# -*- coding: utf8 -*-

import matplotlib.pyplot as plt
import numpy as np

pfiles=['electron','positron','proton','neutron','photon','muon']

atot=2.0*np.pi*280.0**2.0
ptype=3
sb3=np.loadtxt('{0}_flux/{0}_fsb3.out'.format(pfiles[ptype]))
name_out='{}-flux-deff.dat'.format(pfiles[ptype])
Evt_sel=np.logical_and(np.mod(sb3[1,:],2)==1.0,sb3[6,:]==0.0)
sb3=sb3[:,Evt_sel]
Nevt=sb3[0,:]
hpat=sb3[1,:]
nbar=sb3[2,:]
einp=sb3[3,:]
emax=sb3[4,:]
edep=sb3[5,:]
htop=sb3[6,:]
print(np.shape(einp))
table=True
if table==True:
  darray=np.zeros([3,3000000])
  j=0
  for k in range(0,8):
    name='{0}_flux/scicrt1_nt_table_t{1}.csv'.format(pfiles[ptype],k)
    f=open(name,'r')
    for line in f:
      if len(line.split(','))==5.0:
        m=np.fromstring(line,count=5,sep=',')
        darray[0,j]=m[0]
        darray[1,j]=m[2]
        darray[2,j]=np.degrees(m[4])
        j+=1
    f.close()

# para n,p,m de 10MeV a 100GeV
# para ph,pos,e low 10MeV a 1GeV, high 100MeV a 100GeV
ebins=np.ravel(np.outer(10**np.arange(1,4),np.arange(1,10,0.5)))
darray=darray[:,darray[0,:]!=0]
Evt_sel=np.isin(darray[0,:],Nevt)
data_nsel=darray[1,:]
ehist_in,b=np.histogram(data_nsel,bins=ebins)
ehist_out,b=np.histogram(einp,bins=ebins)
dEh=ehist_out/ehist_in
# archivos n,p,m de 0 - 28
#archivos ph,pos,e - low 0 -9, high 0 - 19

fig,ax=plt.subplots(nrows=1,ncols=1,sharex=False,sharey=False)
ax.semilogx(ebins[:-1],dEh,ds='steps')
plt.show()
