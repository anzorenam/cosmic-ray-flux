#!/usr/bin/env python3.7
# -*- coding: utf8 -*-

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integral
import scipy.stats as stats
import seaborn as sns

sns.set(rc={"figure.figsize":(8,4)})
sns.set_context('paper',font_scale=1.5,rc={'lines.linewidth':1.5})
sns.set_style('ticks')
plt.rc('text',usetex=True)
plt.rc('text.latex',preamble=r'\usepackage[utf8]{inputenc} \usepackage[T1]{fontenc} \usepackage[spanish]{babel} \usepackage[scaled]{helvet} \renewcommand\familydefault{\sfdefault} \usepackage{amsmath,amsfonts,amssymb} \usepackage{siunitx}')

pfiles=['electron','positron','proton','neutron','photon','muon']
atot=2.0*np.pi*280**2.0
wilson=True
FdE=np.zeros(28)
dN=0
c0=sns.cubehelix_palette(8,rot=-.4,reverse=True)
fig,ax=plt.subplots(nrows=1,ncols=1,sharex=False,sharey=False)
for pt in range(0,6):
  dI=np.fliplr(np.loadtxt('{0}-spectrum.dat'.format(pfiles[pt])))
  dEff=np.loadtxt('{0}-deff.dat'.format(pfiles[pt]))
  p,x,n=dEff[0,:],dEff[0,:]*dEff[1,:],dEff[1,:]
  abins=np.radians(np.arange(0,90))
  if wilson==False:
    pH,pL=stats.beta.ppf(1-0.025,x+1,n-x),stats.beta.ppf(0.025,x,n-x+1)
  else:
    pL=(2.0*x+1.96**2.0-(1.96*np.sqrt(1.96**2.0-(1.0/n)+4.0*x*(1.0-p)+(4.0*p-2.0))+1.0))/(2.0*(n+1.96**2.0))
    pH=(2.0*x+1.96**2.0+(1.96*np.sqrt(1.96**2.0-(1.0/n)+4.0*x*(1.0-p)-(4.0*p-2.0))+1.0))/(2.0*(n+1.96**2.0))
  ebins=np.ravel(np.outer(10**np.arange(1,5),np.arange(1,10)))
  for j in range(0,28):
    FdE[j]=integral.trapz(np.sin(abins)*dI[j,:],abins)
  Fint=integral.trapz(p*FdE[0:28],ebins[0:28])
  if pt==5:
    Fint=2.1052057266710884*Fint
  ax.errorbar(ebins[0:28],p,yerr=(p-pL,pH-p),color=c0[pt-2])
  ax.set_xscale('log')
  ax.set_yscale('log')
  dN+=atot*Fint
print(dN*3600)
plt.xlabel(r'Energy $[\si{\mega\electronvolt}]$',x=0.9,horizontalalignment='right')
plt.ylabel(r'Detection efficiency')
plt.xlim(5e0,1e4)
plt.ylim(1e-5,1e0)
plt.tight_layout(pad=1.0)
plt.show()
