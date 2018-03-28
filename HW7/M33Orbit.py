#Jenny Calahan
#Astr 400B
#Hw 7 predicting the future of M33 as it orbits M31

import numpy as np
import astropy.units as u
from ReadFile import Read

class M33AnalyticOrbit:
    
    def __init__(self,filename):
        
        #Define universal constants
        self.G= 4.498768e-6
        
        #Define x, y, z postion of M33 as compared to M31 (taken from HW4)
        self.x = -98. # -98.56  
        self.y = -120. #-119.98 
        self.z = -127. #-127.76 
        
        #Define vx, vy, vz of M33 as compared to M31 (taken from HW4)
        self.vx = -28. #-28.41 
        self.vy = 174. #173.92
        self.vz = 92. #93.22
        
        #Disk, Bulge, and Halo data of M31 taken from HW3 all in kpc
        self.rd = 5. 
        self.Mdisk = 0.12e12
        self.rbulge = 1.
        self.Mbulge = 0.019e12
        self.rhalo = 62.            #computed in HW5
        self.Mhalo = 1.921e12
        
    def HernquistAccel(self,M,r_a,x,y,z,i):  
        #Inputs: M = total halo/bulge mass
        #        r_a = scale length
        #        x,y,z = position components 
        #        i = the component of acceleration the user is looking for (x, y, or z)
        
        if i=='z':
            c=z
        if i=='y':
            c=y
        if i=='x':
            c=x
        
        r=np.sqrt(x**2+y**2+z**2)         #find the position vector
        a=-(self.G*M*c)/(r*(r_a+r)**2)    #find the hernquist acceleration
        
        #Output: hernquist acceleration
        return a
    
    def MiyamotoNagaiAccel(self,M,rd,x,y,z,i):
        #Inputs: M = total disk mass
        #        rd = disk scale length
        #        x,y,z = position components 
        #        i = the component of acceleration the user is looking for (x, y, or z)        
        
        zd=rd/5.               #define disk scale height
        B=rd+np.sqrt(z**2+zd**2)   
        R=np.sqrt(x**2+y**2)       
        
        if i == 'z':
            a=-(self.G*M*B*z)/(((R**2+B**2)**1.5)*(np.sqrt(z**2+zd**2)))
        if i == 'y':
            a=-(self.G*M*y)/((R**2+B**2)**1.5)
        if i == 'x':
            a=-(self.G*M*x)/((R**2+B**2)**1.5)
            
        #Output: acceleration as calculated by Miyamoto-Nagai profile and component
        return a
    
    def M31Accel(self,x,y,z,i):
        #Input x,y,z = position
        #      i = direction (x, y, or z)
        
        #calculate bulge acceleration in i direction
        a_b=self.HernquistAccel(self.Mbulge,self.rbulge,x,y,z,i)
        #calculate halo acceleration in i direction
        a_h=self.HernquistAccel(self.Mhalo,self.rhalo,x,y,z,i)
        #calculate disk acceleration in i direction
        a_d=self.MiyamotoNagaiAccel(self.Mdisk,self.rd,x,y,z,i)
        
        tot=a_b+a_h+a_d
        #Output: the sum of the accelerations and direction i
        return tot
        
    
    def LeapFrog(self,dt,x,y,z,vx,vy,vz):
        #Inputs: dt = time interval for integration
        #        x,y,z = starting position vector of M33's COM
        #        vx,vy,vz = starting velocity vector of M33's COM
        
        #Define new x, y, and z postions after half a time step
        xh=x+vx*dt/2.
        yh=y+vy*dt/2.
        zh=z+vz*dt/2.
        
        #Find the x, y, and z acceleration of the galaxy
        M31Ax=self.M31Accel(xh,yh,zh,'x')
        M31Ay=self.M31Accel(xh,yh,zh,'y')
        M31Az=self.M31Accel(xh,yh,zh,'z')
        
        #calulate velcoity at a full time step
        vxn=vx+M31Ax*dt
        vyn=vy+M31Ay*dt
        vzn=vz+M31Az*dt
        
        #calculate final position at a full time step
        xn=x+.5*(vx+vxn)*dt
        yn=y+.5*(vy+vyn)*dt
        zn=y+.5*(vz+vzn)*dt
        
        #Output: an array containing position and velocity components
        return [xn,yn,zn,vxn,vyn,vzn]
    
    def OrbitIntegrator(self,t0,tf,dt):
        #Inputs: t0 = start time
        #        dt = time interval
        #        tf = final time
        
        #Define starting positions:
        x=self.x
        y=self.y
        z=self.z
        vx=self.vx
        vy=self.vy
        vz=self.vz
        
        #Define empty 2D array to hold final position and velocity information
        steps=int((tf-t0)/dt)+1
        Orbit = np.zeros((steps,7)) 
        
        Orbit[0,1],Orbit[0,2],Orbit[0,3]=x,y,z
        Orbit[0,4],Orbit[0,5],Orbit[0,6]=vx,vy,vz
        Orbit[0,0]=t0
        
        np.savetxt('M33Ana', Orbit, header='t   x    y    z   vx   vy   vz', comments='# ',
                fmt=['%.2f', '%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])
        t=t0+dt
        n=1
        while t <= tf:
            #call Leap Frog to integrate over one step
            values = self.LeapFrog(dt,x,y,z,vx,vy,vz)

            #Assign calculated position vectors to empty arrays
            Orbit[n,1],Orbit[n,2],Orbit[n,3]=values[0],values[1],values[2]
            
            #Assign calculated velocity vectors to empty arrays
            Orbit[n,4],Orbit[n,5],Orbit[n,6]=values[3],values[4],values[5]
            
            #Assign time to Orbit
            Orbit[n,0]=t
            
            #Set up new variables
            n=1+n
            t=t+dt
            x,y,z=values[0],values[1],values[2]
            vx,vy,vz=values[3],values[4],values[5]
        
        #Save results into text file
        np.savetxt('M33Ana', Orbit, header='t   x    y    z   vx   vy   vz', comments='# ',
                fmt=['%.2f', '%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])
        
        return n
        
        
        
        
        
        
        
        
        
        
        
        
        
        