#Jenny Calahan
#Astr 400B
#Febuary 18, 2018

#Store the separation and relative velocities of the simulated galaxies during the entire simulation

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from ReadFile import Read
from CenterofMass import *

def OrbitCOM(galaxy,start,end,n):
    #Inputs: Galaxy name, start and end times and n # of space inbetween start and end
    THEfilename = "Orbit_galaxyname.txt"
    Orbit = np.zeros((int(end/n)+1,7))
    
    #For M33
    #delta = 0.5
    #VolDec = 4

    #For MW
    #delta = 1
    #VolDec = 2

    #For M31
    delta = 1
    VolDec = 2

    for i in np.arange(start,end,n):    
        #add a string to the filenumber to the value '000'
        ilbl = '000' + str(i)
        #remove all but the last three digits
        ilbl = ilbl[-3:]
        filename="%s_"%(galaxy) + ilbl + '.txt'

        #Create a center of mass object
        COM = CenterOfMass('VLowRes/'+filename, 2)

        #Store time in the first column in Gyr
        Orbit[int(i/n),0]=COM.time.value/1000

        #Store Position Vector of COM
        pos = COM.COM_P(delta,VolDec)
        Orbit[int(i/n),1]=pos[0]
        Orbit[int(i/n),2]=pos[1]
        Orbit[int(i/n),3]=pos[2]
        
        #Store Velocity Vector of COM
        vel = COM.COM_V(delta,VolDec)
        Orbit[int(i/n),4]=vel[0].value
        Orbit[int(i/n),5]=vel[1].value
        Orbit[int(i/n),6]=vel[2].value
        
        print(i)

    #create file for Orbit information to go in
    fileout='M31_lifetime.txt'
    np.savetxt(fileout, Orbit, header='t x y z vx vy vz', comments='#',fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])
    #Output: File galaxy_lifetime.txt and Orbit, 2D array that goes into the file
    return Orbit

    
        

        














