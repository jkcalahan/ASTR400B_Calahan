#Jenny Calahan
#Astr 400B HW5
#Feb 13th, 2018

#import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
from astropy import constants as const
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from ReadFile import *
from CenterofMass import *

class MassProfile:

    def __init__(self, galaxy, snap):

        #add a string to the filenumber to the value '000'
        ilbl = '000' + str(snap)
        #remove all but the last three digits
        ilbl = ilbl[-3:]
        self.filename="%s_"%(galaxy) + ilbl + '.txt'

        # read in the file and particle type
        self.time, self.total, self.data = Read(self.filename)
        
        #create an array to store indexes of particles of desired Ptype
        #self.index = np.where(self.data['type'] == ptype)
        self.index = range(len(self.data))

        # store the mass, positions, velocities of only the particles of the given type
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]
        self.haloind = np.where(self.data['type'] == 1)
        self.mhalo = self.data['m'][self.haloind]
        self.gname = galaxy

    def MassEnclosed(self,ptype,radii):     #Inputs: particle type and an array of radii
        #Find the center of mass position of a galaxy using particles of ptype
        if self.gname == 'M33' and ptype == 3:
            masses=np.zeros(len(radii))
        else:
            COM = CenterOfMass(self.filename, ptype)
            COMpos = COM.COM_P(1)

            #create arrays of positions and mass of objects of ptype
            ind = np.where(self.data['type'] == ptype)
            newx=self.x[ind]
            newy=self.y[ind]
            newz=self.z[ind]
            m = self.data['m'][ind]

            #Create an array of distances that objects are from the center of mass
            dist=np.sqrt(((newx-COMpos[0])**2+(newy-COMpos[1])**2+(newz-COMpos[2])**2))*u.kpc

            masses=[]                                      #an array to hold the mass enclosed at certain radii
            for r in radii:
                masstemp=0
                for dind in range(len(dist)):
                    if dist[dind] < r*u.kpc:
                        masstemp=masstemp+m[dind]*1e10         #add up all the masses of objs less than r
                masses.append(masstemp)                    #when we go though all objs, store final mass in array masses

        #Outputs: an array of masses at different radii  
        return masses

    def MassEnclosedTotal(self,radii):    #Input: an array of radii that we want to know the tot mass enclosed 
        #Use MassEnclosed to find disk, bulge, and halo mass. 
        #Only excpetion is if galaxy is M33, then bulge mass is set to zero
        bulge_mass=self.MassEnclosed(3,radii)
        disk_mass=self.MassEnclosed(2,radii)
        halo_mass=self.MassEnclosed(1,radii)

        #Find the total mass of the galaxy at each radii
        totmass=[]
        for j in range(len(radii)):
            totmass.append((bulge_mass[j]+disk_mass[j]+halo_mass[j]).value)
        #Output: array of masses of all particles at different radii
        return totmass

    def HernquistMass(self,radii,a,Mhalo):    #Inputs: radius (arr), scaling factor(float), mass of the halo (arr)
        mass=[]
        #Sooooo, I don't use the user's Mhalo, I calculate the total halo mass internally: self.mhalo, see __init__
        for R in radii:
            mass.append((sum(self.mhalo)*1e10*R**2)/(a+R)**2)
 
        return mass                       #Output: Halo mass 

    def CircularVelocity(self,ptype,radii):  #Inputs: particle type and array of radii
        #find the mass enclosed at certain radii
        if self.gname =='M33' and ptype == 3:
            cirvel=np.zeros(len(radii))
        else:
            masses=self.MassEnclosed(ptype,radii)
        
            #Calculate the circular velocity by using vc=sqrt(GM/r)
            cirvel=[]
            G=const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
            for j in range(len(radii)):
                cirvel.append((np.sqrt(G*masses[j]/radii[j]/u.kpc)).value)

        #Output: array of circular velocities at each radii 
        return cirvel

    def CircularVelocityTotal(self,radii):
        #Find the circular velocity of each component at each radii
        totcircvel=[]
        G=const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        #Find the mass enclosed for each component then add each component at the same radii to find circ vel
        if self.gname =='M33':
            halom=self.MassEnclosed(1,radii)
            diskm=self.MassEnclosed(2,radii)
            for j in range(len(halom)):totcircvel.append(np.sqrt(G*(halom[j]+diskm[j])/radii[j]/u.kpc).value)
        else:
            halom=self.MassEnclosed(1,radii)
            diskm=self.MassEnclosed(2,radii)
            bulgem=self.MassEnclosed(3,radii)
            for j in range(len(halom)):totcircvel.append(np.sqrt(G*(halom[j]+diskm[j]+bulgem[j])/radii[j]/u.kpc).value)

        #Output: An array of total circular velocity at each radii
        return totcircvel

    def HernquistVCirc(self,radii,a,Mhalo):    #Inputs: radius, scaling factor, mass of the halo
        mass=self.HernquistMass(radii,a,Mhalo)
        #Output: circular velocity using mass as caluclated by the Hernquist mass function
        G=const.G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        vcirc=[]
        for i in range(len(mass)): vcirc.append(np.sqrt(G*mass[i]/radii[i]/u.kpc).value)
        return vcirc



        
        

        











