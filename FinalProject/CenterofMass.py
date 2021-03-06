# Homework 4 
# Center of Mass Position and Velocity

#Edits made by Jenny Calahan
#Feb 04, 2018

#Edit made on Feb 17th: in Com_P changed RMAX to be divided by user input VolDec instead of 2.0

# import modules
import numpy as np
import astropy.units as u
from ReadFile import Read


class CenterOfMass:
   
    def __init__(self, filename, ptype):
        # read in the file 
        self.time, self.total, self.data = Read(filename)
        #print(self.time)
            
        #create an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)
    
        # store the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]

    
    
    def COMdefine(self,i,j,k,m):
    # Function to compute the center of mass position or velocity generically
    # input: array (a,b,c) of positions or velocities and the mass
    # returns: 3 floats  (the center of mass coordinates)

        # xcomponent Center of mass  
        Icom = np.sum(i*m)/np.sum(m)
        # ycomponent Center of mass
        Jcom = np.sum(j*m)/np.sum(m)
        # zcomponent
        Kcom = np.sum(k*m)/np.sum(m)

        return Icom, Jcom, Kcom

        
    def COM_P(self, delta,VolDec):
    # Function to specifically return the center of mass position and velocity
    # input:    
    #        particle type (1,2,3)
    #        delta (tolerance)
    # returns: One vector, with rows indicating:
    # time in Gyr 
    # 3D coordinates of the center of mass position (kpc), 
    # 3D velocity vector of the center of mass (km/s)

       

        # Center of Mass Position
        ###########################
        
        # Try a first guess at the COM position by calling COMdefine
        XCOM, YCOM, ZCOM = self.COMdefine(self.x,self.y,self.z,self.m)
        # compute the magnitude of the COM position vector. 
        RCOM = np.sqrt(XCOM**2 + YCOM**2 + ZCOM**2)
        # print('init R', RCOM)

    
        # iterative process to determine the center of mass
    
        # change reference frame to COM frame
        # compute the difference between particle coordinates
        # and the first guess at COM position  
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        RNEW = np.sqrt(xNew**2.0 + yNew**2.0 +zNew**2.0)
    
        # find the max 3D distance of all particles from the guessed COM
        # will re-start at half that radius (reduced radius)
        RMAX = max(RNEW)/VolDec

        # pick an initial estimate for the change in COM position 
        # between the first guess above and the new one computed from half that volume.
        CHANGE = 1000.0
    
        # start iterative process to determine center of mass position
        # delta is the tolerance for the difference in the old COM and the new one.
        while (CHANGE > delta):
     
            # select all particles within the reduced radius (starting from original x,y,z, m)
            index2 = np.where(RNEW < RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]
        
            # Refined COM position:
            # compute the center of mass position using 
            # the particles in the reduced radius
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2,y2,z2,m2)
            # compute the new 3D COM position
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)

            # determine the difference between the previous center of mass position 
            # and the new one. 
            CHANGE = np.abs(RCOM - RCOM2)
            # check this
            # print ("DIFF", diff)
        
            # Before loop continues, reset : RMAX, particle separations and COM
        
            # reduce the volume by a factor of 2 again
            RMAX = RMAX/VolDec
            # check this.
            #print ("maxR", maxR)
        
        
            # Change the frame of reference to the newly computed COM.
            # subtract the new COM
            xNew = self.x - XCOM2
            yNew = self.y - YCOM2
            zNew = self.z - ZCOM2
            RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)
        
      
        
            # set the center of mass positions to the refined values
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2
     
            # create a vector to store the COM position 
            # set the correct units usint astropy   
            # round all values
            #COMP = [np.round((XCOM)*u.kpc), np.round((YCOM)*u.kpc), np.round((ZCOM)*u.kpc)]
            COMP = [np.round((XCOM)), np.round((YCOM)), np.round((ZCOM))]
    
        # return the COM positon vector
        return COMP
        
    
    def COM_V(self, delta, VolDec):

        # the max distance from the center that we will use to determine the center of mass velocity
        RVMAX = 15.0*u.kpc
    
        # Determine the center of mass position
        COMP = self.COM_P(delta,VolDec)
              
        # determine the position of all particles relative to the center of mass position

        xV = self.x[:]*u.kpc - COMP[0]*u.kpc
        yV = self.y[:]*u.kpc - COMP[1]*u.kpc
        zV = self.z[:]*u.kpc - COMP[2]*u.kpc
        RV = np.sqrt(xV**2 + yV**2 + zV**2)
    
        # determine the index for those particles within the max radius
        indexV = np.where(RV < RVMAX)
    
        # determine the velocity and mass of those particles within the mas radius
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew = self.m[indexV]
    
        # compute the center of mass velocity using those particles
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew,vynew,vznew, mnew)

        # create a vector to store the COM velocity
        # set the correct units usint astropy   
        # round all values
        COMV = [np.round((VXCOM)*u.km/u.s), np.round((VYCOM)*u.km/u.s), np.round((VZCOM)*u.km/u.s)]
    
        # return the COM vector
        return COMV
