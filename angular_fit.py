#!/usr/bin/env python3.7
# -*- coding: utf8 -*-

import matplotlib as mat
import matplotlib.pyplot as plt
import numpy as np
import ROOT
import root_numpy as rnp

pflux=1e-2*np.array([2.1873682915465129,1.7103582799966752,0.96944022115073025,
                     5.8755789311196883,14.578720064501247,3.48302386565864360])

pfiles=['electron','positron','proton','neutron','photon','muon']

fit=False
flux=True
if fit==True:
  theta=np.zeros([6,100000])
  zdist=np.zeros([6,90])
  abins=np.arange(0,90)
  abins_rad=np.radians(abins)
  for k in range(5,6):
    print(pfiles[k])
    name='{0}_sn.out'.format(pfiles[k])
    f=open(name,'r')
    j=0
    for line in f:
      theta[k,j]=np.float(line.split()[2])
      j+=1
    theta[k,:]=np.degrees(np.arccos(theta[k,:]))
    angraph=ROOT.TGraph()
    ang_hist,bin=np.histogram(theta[k,:],bins=abins,density=True)
    angT=np.transpose(np.array([abins_rad[:-1],ang_hist]))
    rnp.fill_graph(angraph,angT)
    if pfiles[k]=='muon':
      afun=ROOT.TF1('angular','[0]*sin(x)*TMath::Power(cos(x),[1])',0,np.pi/2.0)
    else:
      afun=ROOT.TF1('angular','[0]*sin(x)*(1.0+[1]*TMath::Power(cos(x),[2]))',0,np.pi/2.0)
    afun.SetParLimits(2,0.5,4.0)
    angraph.Fit(afun,'R')
    p=afun.GetParameters()
    e=afun.GetParError(0)
    if pfiles[k]=='muon':
      zdist[k,:]=p[0]*np.sin(abins_rad)*(np.power(np.cos(abins_rad),p[1]))
      print(p[0],p[1])
    else:
      zdist[k,:]=p[0]*np.sin(abins_rad)*(1.0+p[1]*np.power(np.cos(abins_rad),p[2]))
      print(p[0],p[1],p[2])
    fig,ax=plt.subplots(nrows=1,ncols=1,sharex=False,sharey=False)
    ax.plot(abins,zdist[k,:])
    ax.hist(theta[k,:],bins=abins,density=True,alpha=0.5)
  plt.show()

# muon 0.05570501105665601+/-0.000599031, 2.2562596147377127+/-0.0331894
# neutron 0.009324194806673637+/-0.000123081, 2.8084413116635183+/-0.0436179, 2.1559001756314506+/-0.000991363
# gamma 0.0023213730909956826+/-0.000113165, 26.526169391742343+/-1.29808, 3.0940580647758016+/-0.0313797
# proton 0.003977277048029226+/-0.000183094, 14.89280789753243+/-0.695653, 3.4344517415983535+/-0.00105636
# electron 0.003869391237031857+/-0.000142957, 13.838381804716969+/-0.522228, 2.9744341245359966+/-0.00111054
# positron 0.001985670696803895+/-8.48726e-05, 31.23550599413086+/-1.3326, 3.088193923966495+/-0.00104052
