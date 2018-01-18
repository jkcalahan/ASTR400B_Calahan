#Jenny Calahan
#Jan 16, 2018
#Astr 400B, Hw 2

import numpy as np
import astropy.units as u 
from ReadFile import *

filename='MW_000.txt'

def ParticleInfo(j,filename):           #Returns distance, velocity, and mass of particle at position 'j'
    t,N,data=Read('MW_000.txt')
    index=np.where(data['type']==2)      #Only looking at disk stars
    x=data['x'][index]
    y=data['y'][index]
    z=data['z'][index]
    vx=data['vx'][index]
    vy=data['vy'][index]
    vz=data['vz'][index]
    m=data['m'][index]
    dist=np.sqrt((x[j]**2)+(y[j]**2)+(z[j]**2))       #calc 3D distance
    vel=np.sqrt((vx[j]**2)+(vy[j]**2)+(vz[j]**2))     #calc 3D velocity
    mass=m[j]*1e10
    return dist,vel,mass

j=101
d,v,ma=ParticleInfo(j,'MW_000.txt')     #find the distance(kpc), velocity(km/s), and mass(M_sun) 
                                                 #of the 100th particle in 'MW_000.txt'
dist=np.around(d,3)
vel=np.around(v,3)
mass=np.around(ma,3)

distly=u.kpc.to(u.lightyear,value=dist)    #convert kpc to ly

#Print final Results
print('3D Distance of object %d: %.3f ly'%(j,distly))
print('3D velocity of object %d: %.3f km/s'%(j,vel))
print('Mass of object %d: %.3f solar masses'%(j,mass))


    