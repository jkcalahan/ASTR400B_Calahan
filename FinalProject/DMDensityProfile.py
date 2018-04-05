#Jenny Calahan
#Determining the density profiles of dark matter particles in M33


import numpy as np
from ReadFile import Read
from CenterofMass import CenterOfMass


#Goals -- as of April 4th: create a density profile at snap 0000, determine its slope
#                          find a way to repeat that for different snaps and plot that slope value at each point
#      Functions to create:   DensProf --> outputs array of mass values at different radii, and array of radii
#                             Slope --> ouputs the slope values of density profile 

#read in M31 dark matter particles with in 20 kpc

def ReadMass(snap,limit):
    #Inputs: snap -- time step we want to study M33 at
    #        limit -- distance in kpc from the COM that defines the core of M33
    #
    #Outputs: m -- float, total dark matter mass within limit
    #---------------------------------------
    
    # Determine Filename
    # add a string of the filenumber to the value "000"
    ilbl = '000' + str(snap)
    # remove all but the last 3 digits
    ilbl = ilbl[-3:]
    # create filenames
    filename='%s_'%('M33') + ilbl + '.txt'
    path='/home/jcalahan/VLowRes/'+filename

    # read in the file 
    time, total, data = Read(path)

    #create an array to store indexes of dark matter particles
    index = np.where(data['type'] == 1)

    # store the mass and positions of only the particles of the given type
    m = data['m'][index]
    x = data['x'][index]
    y = data['y'][index]
    z = data['z'][index]
    
    #Find COM postion of dark matter particles
    COM = CenterOfMass(path,1)
    COMP = COM.COM_P(.2,4.0)           #<--using values from HW6
    
    # store distances of particles from COM
    dist=np.sqrt((x-COMP[0])**2+(y-COMP[1])**2+(z-COMP[2])**2)
            
    # initilize tot mass of dark matter particles        
    mass = 0
    
    # 1 particle of dark matter mass (1e10)
    dm_mass = 0.00394985
    
    #if a particle is within the limit, add a dark matter mass to total mass
    for d in dist:
        if d <= limit:
            mass = mass + dm_mass
                  
    return mass*(1e10)



    
    