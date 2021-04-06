#!/usr/bin/env python3.7
# -*- coding: utf8 -*-

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import time
import datetime

def show_track(yside,xside):
  fig,ax=plt.subplots(nrows=2,ncols=1,sharex=False,sharey=False)
  sns.heatmap(yside,vmin=0,vmax=7.0,cmap='viridis',rasterized=True,
              xticklabels=False,yticklabels=False,ax=ax[0])
  plt.title('X-side')
  sns.heatmap(xside,vmin=0,vmax=7.0,cmap='viridis',rasterized=True,
              xticklabels=False,yticklabels=False,ax=ax[1])
  plt.title('Y-side')
  plt.show()

def read_data(data_file,counter):
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
    side=np.uint8(k[4])
    col=np.uint8(k[5]+8*k[6])
    reg=np.uint8(k[3])
    if n0!=nevent:
      n0=nevent
      img=img[:,:,index]
      #yside=img[0,:,:]
      #xside=img[1,:,:]
      yside=np.flipud(img[0,:,:])
      xside=np.flipud(img[1,:,:])
      trg_evt=np.any(yside>0) or np.any(xside>0)
      trg_top=np.any(yside[0,:]>=Emip) or np.any(xside[0,:]>=Emip)
      trg_sb1=np.any(yside[1:8,3:115]>=Ethr) and np.any(xside[1:8,3:115]>=Ethr)
      trg_sb2=np.any(yside[8:16,3:115]>=Ethr) and np.any(xside[8:16,3:115]>=Ethr)
      trg_sb3=np.any(yside[16:24,3:115]>=Ethr) and np.any(xside[16:24,3:115]>=Ethr)
      if trg_evt:
        if trg_sb3:
          sb3_y=yside[16:24,3:115]
          sb3_x=xside[16:24,3:115]
          print('hit pattern: {0}'.format(4.0*trg_sb1+2.0*trg_sb2+1.0*trg_sb3))
          print('anti-coin: {0}'.format(1.0*trg_top))
          show_track(yside,xside)
        counter+=1
      img=np.zeros((2,24,118))
      img[side,col,reg]=k[2]*(k[2]>Emin)
    else:
      img[side,col,reg]=k[2]*(k[2]>Emin)
  f.close()
  return counter

t0=time.time()
workers=np.arange(8)
jevt=0
for w in workers:
  name='scicrt1_nt_edep_t{0}.csv'.format(w)
  jevt=read_data(name,jevt)
t1=time.time()
