#!/bin/python

from os import system,getenv
import subprocess
import matplotlib.pyplot as plt
import numpy as np

mode=1
temp=20
rain=5
pop=100
nstep=500

tempv=np.linspace(17,39,num=80)

vars=["eir","cspr","PRd","cases","immunity","waterfrac"]
varsf=[365,100,100,1000,100,100]
varsu=[" (bites per year)"," (%)"," (%)"," (per thousand)"," (%)"," (%)"]


vals=np.zeros((len(vars),len(tempv)))

for itemp,temp in enumerate(tempv) :
    system("{0}/scripts/vectri_driver {1} {2} {3} {4} {5}".format(getenv("VECTRI"),mode ,temp,rain,pop,nstep))
    system("cdo -timmean output/vectri.nc output/tm.nc")

    for ivar,var in enumerate(vars) : 
        val=subprocess.check_output("ncdump -v,{0} output/tm.nc | tail -n 2 | head -n 1 | awk '{{print $1;}}'".format(var), shell=True)
        print (var,val)
        vals[ivar,itemp]=float(val)*varsf[ivar]

plt.figure(1)
for ivar,var in enumerate(vars): 
    plt.subplot(len(vars)/2,2,ivar)
    plt.plot(tempv,vals[ivar,:])
    plt.xlabel('Temperature (C)')
    plt.ylabel(vars[ivar]+varsu[ivar])
    plt.tight_layout()


plt.savefig('eir_equil.pdf',dpi=300)
plt.clf()


