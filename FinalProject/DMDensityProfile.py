#Jenny Calahan
#Determining the density profiles of dark matter particles in M33


import numpy as np
from ReadFile import Read

#read in M31 dark matter particles with in 20 kpc

def Read(filename)
    # read in the file 
    time, total, data = Read(filename)

    #create an array to store indexes of particles of desired Ptype
    index = np.where(self.data['type'] == ptype)

    # store the mass, positions, velocities of only the particles of the given type
    m = self.data['m'][index]
    x = self.data['x'][index]
    y = self.data['y'][index]
    z = self.data['z'][index]
    vx = self.data['vx'][index]
    vy = self.data['vy'][index]
    vz = self.data['vz'][index]
    
    