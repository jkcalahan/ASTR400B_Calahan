# Homework 4 
# Center of Mass Position and Velocity

#Edits made by Jenny Calahan
#Feb 04, 2018

# import modules
import numpy as np
import astropy.units as u
from ReadFile import *


class CenterOfMass:
   
    def __init__(self, filename, ptype):
        # read in the file and particle type
        self.time, self.total, self.data = Read(filename)
            
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
      
    def total_mass(self,mass):
        #Note: you can add other keyword arguments into the function, but 'self' must be first        
        return np.sum(mass)*u.Msun*1e10

    def COMdefine(self,i,j,k,mass):                      #Goal: determine the x, y, and z components of the center of mass
        topi,topj,topk=0,0,0                             #Eventually top will equal the sum of x_i*m_i
        for c in range(len(i)):     topi=topi+(i[c]*mass[c])       #summing the product of x pos and mass
        for c in range(len(j)):     topj=topj+(j[c]*mass[c])       #summing the product of y pos and mass
        for c in range(len(k)):     topk=topk+(k[c]*mass[c])       #summing the product of z pos and mass
        topi=topi*u.Msun*1e10                            #put sum of masses and positions in same units of tot. mass
        topj=topj*u.Msun*1e10
        topk=topk*u.Msun*1e10
        bottom=self.total_mass(mass)                     #bottom of the fraction will be total mass
        return topi/bottom, topj/bottom, topk/bottom     #x, y, and z components of the center of mass of any galaxy

    def COM_P(self,delta):    #Input: Delta value to determine COM
        xtemp,ytemp,ztemp,mtemp=self.x,self.y,self.z,self.m     #start with an inital array of objects' x,y, and z coordinates
        RCOMdiff=delta+10                                       #define an arbitrary RCOM that is larger than delta
        #continue to recalculated RCOM until  the difference is less than delta
        while RCOMdiff > delta:
            x1,y1,z1=self.COMdefine(xtemp,ytemp,ztemp,mtemp)  #pull the initial COM coordinates
            xnew,ynew,znew=[],[],[]                        #define new arrays for objects' x, y, and z coordinates in COM frame
            RNEW=[]                                        #an array of distances of the objects in COM frame
            RCOM=np.sqrt(x1**2+y1**2+z1**2)                #distance of the center of mass
            for j in xtemp: xnew.append(j-x1)              #subract XCOM from each x values and put in xnew as defined above
            for j in ytemp: ynew.append(j-y1)              #^^ same but in y, and then z
            for j in ztemp: znew.append(j-z1)

            #xtemp ytemp and ztemp contain x, y, and z coordinates in COM frame 1
            for k in range(len(xnew)): RNEW.append(np.sqrt(xnew[k]**2+ynew[k]**2+znew[k]**2))  #array of distances in COM frame
            RMAX=max(RNEW)/2.                            #the max 3D seperation of an obj from COM divided by 2

            #Now, we're going to caluclate a new COM using only objs less than RMAX away from the initial COM
            xtemp2,ytemp2,ztemp2,mtemp2=[],[],[],[]
            for l in range(len(xnew)):                   #make new arrays, xtemp2 etc. that will be objs within RMAX
                if np.sqrt(xnew[l]**2+ynew[l]**2+znew[l]**2) < RMAX:
                    xtemp2.append(xtemp[l])               #if an object is less than RMAX away, it will go into our calculate
                    ytemp2.append(ytemp[l])               #of a new COM
                    ztemp2.append(ztemp[l])               #Generally speaking 'temp' will be in galactic frame 'new in COM frame
                    mtemp2.append(mtemp[l])

            #A new COM using only objs less than RCOM away from the inital COM:
            x2,y2,z2=self.COMdefine(xtemp2,ytemp2,ztemp2,mtemp2)    
            RCOM2=np.sqrt(x2**2+y2**2+z2**2)      #New distance to the new COM
            RCOMdiff=np.absolute(RCOM-RCOM2)      #Difference between initial COM and COM of objs > RMAX
            xtemp,ytemp,ztemp,mtemp=xtemp2,ytemp2,ztemp2,mtemp2   #in case we need to run RCOMdiff > delta, we have new xtemp

        #Output: Once the difference between the two COMs is less than delta, print out the most recent calc of COM
        return x2,y2,z2           
    
    def COM_V(self,xcom,ycom,zcom):   #Inputs: center of mass x, y, and z coordinates
        dist,xvel,yvel,zvel,marray=[],[],[],[],[]
        for i in range(len(self.x)):            #Make an array of distances of objs in the COMframe
            dist.append(np.sqrt((self.x[i]-xcom)**2+(self.y[i]-ycom)**2+(self.z[i]-zcom)**2))
        for j in range(len(dist)):
            if dist[j] < 15:          #if an obj is located less than 15 kpc, then they're velocities are added to xvel...
                xvel.append(self.vx[j])
                yvel.append(self.vy[j])
                zvel.append(self.vz[j])
                marray.append(self.m[j])
        vxcom,vycom,vzcom=self.COMdefine(xvel,yvel,zvel,marray)   #find the com velocity with only objs <15 kpc
        return vxcom,vycom,vzcom       #Outputs: the x, y, and z components of the velocity of the center of mass

# EXAMPLE OF USING A CLASS
##########################

filename, ptype='M31_000.txt',2
t,N,data=Read(filename)                     #t == time, N==number of particles, data==all data
index=np.where(data['type']==ptype)
m=data['m'][index]
x=data['x'][index]
y=data['y'][index]
z=data['z'][index]
vx=data['vx'][index]
vy=data['vy'][index]
vz=data['vz'][index]

# Create a Center of mass object for the MW
MWCOM = CenterOfMass('MW_000.txt', 2)
M31COM = CenterOfMass('M31_000.txt',2)
M33COM = CenterOfMass('M33_000.txt',2)
DELTA=.3           #delta is 300 pc

#MWCOMx,MWCOMy,MWCOMz = MWCOM.COM_P(DELTA)
#print('The center of mass of the MW is %.3f x_hat %.3f y_hat and %.3f z_hat kpc'%(MWCOMx,MWCOMy,MWCOMz))
#MWCOMvx,MWCOMvy,MWCOMvz = MWCOM.COM_V(MWCOMx,MWCOMy,MWCOMz)
#print("The center of mass velocity of the MW is: %.3f x_hat %.3f y_hat and %.3f z_hat km/s"%(MWCOMvx,MWCOMvy,MWCOMvz))

#xCOM,yCOM,zCOM=M31COM.COMdefine(x,y,z,m)
#print(xCOM,yCOM,zCOM)

#M31COMx,M31COMy,M31COMz = M31COM.COM_P(DELTA)
#print('The center of mass of M31 is %.3f x_hat %.3f y_hat and %.3f z_hat kpc'%(M31COMx,M31COMy,M31COMz))
#M31COMvx,M31COMvy,M31COMvz = M31COM.COM_V(M31COMx,M31COMy,M31COMz)
#print(M31COMvx,M31COMvy,M31COMvz)

