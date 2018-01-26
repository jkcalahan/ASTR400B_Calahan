#Jenny Calahan
#Astr 400B
#HW 3
#Jan 23rd, 2018

#Import Modules 
import numpy as np
import astropy.units as u
from ReadFile import *

def ComponentMass(filename,ptype,):             #returns the total mass of any desired galaxy component
    t,N,data=Read(filename)                     #t == time, N==number of particles, data==all data
    index=np.where(data['type']==ptype)
    m=data['m'][index]
    mtemp=0
    for i in m:                                 #Add all individual masses of type ptype objs
        mtemp=i+mtemp
    TotMass=np.around(mtemp,3)*1e10             #Mulitply mass by 1e10 to get in terms of solar mass and round by 3 dec places    

    return TotMass                              #return total mass of all ptype objects

#Tot Mass of objs in MW in units of 10^12 sol mass:
MWHalo=ComponentMass('MW_000.txt',1)/(1e12)
MWDisk=ComponentMass('MW_000.txt',2)/(1e12)
MWBulge=ComponentMass('MW_000.txt',3)/(1e12)
MWtot=MWHalo+MWDisk+MWBulge
MWf=(MWDisk+MWBulge)/(MWtot)
print('Milky Way & %.3f & %.3f & %.3f & %.3f & %.3f \\\ '%(MWHalo,MWDisk,MWBulge,MWtot,MWf))  
                                                                              #prints results in LaTeX friendly form

#Tot Mass of objs in M31 in units of 10^12 sol mass:
M31Halo=ComponentMass('M31_000.txt',1)/(1e12)
M31Disk=ComponentMass('M31_000.txt',2)/(1e12)
M31Bulge=ComponentMass('M31_000.txt',3)/(1e12)
M31tot=M31Halo+M31Disk+M31Bulge
M31f=(M31Disk+M31Bulge)/(M31tot)
print('M31 & %.3f & %.3f & %.3f & %.3f & %.3f \\\ '%(M31Halo,M31Disk,M31Bulge,M31tot,M31f))

#Tot Mass of obj in M33 in units of 10^12 sol mass:
M33Halo=ComponentMass('M33_000.txt',1)/(1e12)
M33Disk=ComponentMass('M33_000.txt',2)/(1e12)
M33Bulge=ComponentMass('M33_000.txt',3)/(1e12)
M33tot=M33Halo+M33Disk+M33Bulge
M33f=(M33Disk+M33Bulge)/(M33tot)
print('M33 & %.3f & %.3f & %.3f & %.3f & %.3f \\\ '%(M33Halo,M33Disk,M33Bulge,M33tot,M33f))

#Print out results for Local Group (MW+M31+M33)
LocalHalo=MWHalo+M31Halo+M33Halo
LocalDisk=MWDisk+M31Disk+M33Disk
LocalBulge=MWBulge+M31Bulge+M33Bulge
LocalTot=LocalHalo+LocalDisk+LocalBulge
Localf=(LocalDisk+LocalBulge)/(LocalTot)
print('Local Group & %.3f & %.3f & %.3f & %.3f & %.3f \\\ '%(LocalHalo,LocalDisk,LocalBulge,LocalTot,Localf))








