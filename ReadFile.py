#Jenny Calahan
#Jan 16, 2018
#Astr 400B Homework 2

import numpy as np
import astropy.units as u

def Read(filename):
    file=open(filename,'r')
    line1=file.readline()
    label, value = line1.split()
    time = float(value)*10.0*u.Myr        #Time
    line2=file.readline()
    label2, value2 = line2.split()
    NParticles = float(value2)            #Number of Particles
    file.close()
    data=np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    return time,NParticles,data

t,N,d=Read('MW_000.txt')


    