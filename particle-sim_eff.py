#!/usr/bin/env python3.7
# -*- coding: utf8 -*-

import matplotlib.pyplot as plt
import numpy as np
import argparse
import time
import datetime
import os

def read_data(data_file,earray,counter):
  bi=np.arange(116,-2,-2)
  bd=np.arange(1,118,2)
  index=np.concatenate((bi,bd),axis=0)
  img=np.zeros((2,24,118))
  f=open(data_file,'r')
  for i in range(0,22):
    f.readline()
  line=f.readline()
  k=np.fromstring(line,dtype=np.float,count=7,sep=',')
  n0=k[0]
  e0=k[1]
  Emin=1.635/14
  Ethr=7.0
  Emip=0.4905
  side=np.uint8(k[4])
  col=np.uint8(k[5]+8*k[6])
  reg=np.uint8(k[3])
  img[side,col,reg]=k[2]*(k[2]>Emin)
  for line in f:
    k=np.fromstring(line,dtype=np.float,count=7,sep=',')
    nevent=k[0]
    ein=k[1]
    side=np.uint8(k[4])
    col=np.uint8(k[5]+8*k[6])
    reg=np.uint8(k[3])
    if n0!=nevent:
      img=img[:,:,index]
      yside=np.flipud(img[0,:,:])
      xside=np.flipud(img[1,:,:])
      trg_evt=np.any(yside>0) or np.any(xside>0)
      trg_top=np.any(yside[0,:]>=Emip) or np.any(xside[0,:]>=Emip)
      trg_sb1=np.any(yside[1:8,3:115]>=Ethr) and np.any(xside[1:8,3:115]>=Ethr)
      trg_sb2=np.any(yside[8:16,3:115]>=Ethr) and np.any(xside[8:16,3:115]>=Ethr)
      trg_sb3=np.any(yside[16:24,3:115]>=Ethr) and np.any(xside[16:24,3:115]>=Ethr)
      hit_pattern=4.0*trg_sb1+2.0*trg_sb2+1.0*trg_sb3
      if trg_evt==True:
        earray[0,counter]=n0
        earray[1,counter]=hit_pattern
        sb3_y=yside[16:24,3:115]
        sb3_x=xside[16:24,3:115]
        earray[2,counter]=np.sum(sb3_y>0)+np.sum(sb3_x>0)
        earray[3,counter]=e0
        earray[4,counter]=np.maximum(np.amax(sb3_y),np.amax(sb3_x))
        earray[5,counter]=np.sum(sb3_y)+np.sum(sb3_x)
        earray[6,counter]=1.0*trg_top
        counter+=1
      n0=nevent
      e0=ein
      img=np.zeros((2,24,118))
      img[side,col,reg]=k[2]*(k[2]>Emin)
    else:
      img[side,col,reg]=k[2]*(k[2]>Emin)
  f.close()
  return earray,counter

parser=argparse.ArgumentParser()
parser.add_argument('file_name', help='archivo de entrada',type=str)
parser.add_argument('file_type', help='tipo de archivo (datos)',type=str)
parser.add_argument('file_out', help='salida',type=str)
args=parser.parse_args()
file_name=args.file_name
file_type=args.file_type
file_out=args.file_out

t0=time.time()
dir='/run/media/shirokuma/neutron-data/cosmic-ray_flux/proton_flux'
plot=False
workers=np.arange(8)
jevt=0
nevent=3000000
edata=np.zeros([7,nevent])
for w in workers:
  name='{0}/{1}1_nt_{2}_t{3}.csv'.format(dir,file_name,file_type,w)
  edata,jevt=read_data(name,edata,jevt)

Evt_sel=edata[0,:]!=0
Nevt=edata[0,Evt_sel]
hpat=edata[1,Evt_sel]
nbar=edata[2,Evt_sel]
einp=edata[3,Evt_sel]
emax=edata[4,Evt_sel]
edep=edata[5,Evt_sel]
htop=edata[6,Evt_sel]
ebins=np.arange(1,1e4,1.0)
print('Eventos totales: {0}'.format(jevt))

fdat=open('{0}_fsb3.out'.format(file_out),'w')
np.savetxt(fdat,Nevt,fmt='%1.4f',newline=' ')
fdat.write('\n')
np.savetxt(fdat,hpat,fmt='%1.4f',newline=' ')
fdat.write('\n')
np.savetxt(fdat,nbar,fmt='%1.4f',newline=' ')
fdat.write('\n')
np.savetxt(fdat,einp,fmt='%1.4f',newline=' ')
fdat.write('\n')
np.savetxt(fdat,emax,fmt='%1.4f',newline=' ')
fdat.write('\n')
np.savetxt(fdat,edep,fmt='%1.4f',newline=' ')
fdat.write('\n')
np.savetxt(fdat,htop,fmt='%1.4f',newline=' ')
fdat.close()

if plot==True:
  fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
  n,bins,patch=ax.hist(edep,bins=ebins,histtype='stepfilled',log=True)
  plt.xscale('log')
  plt.tight_layout(pad=1.0)
  plt.savefig('neutron-edep.pdf')

  fig,ax=plt.subplots(nrows=1,ncols=1,sharex=True,sharey=True)
  n,bins,patch=ax.hist(einp,bins=ebins,histtype='stepfilled',log=True)
  plt.xscale('log')
  plt.tight_layout(pad=1.0)
  plt.savefig('neutron-einp.pdf')

t1=time.time()
dtime=datetime.timedelta(seconds=(t1-t0))
print('Tiempo total {0:1.2f} seg.'.format(dtime.total_seconds()))
