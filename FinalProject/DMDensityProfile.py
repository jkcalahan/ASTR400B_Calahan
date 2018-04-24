#Jenny Calahan
#Determining the density profiles of dark matter particles in M33


import numpy as np
from ReadFile import Read
from CenterofMass import CenterOfMass
from scipy.optimize import curve_fit


#Goals -- as of April 21st:   Find the slope of density profiles at various timesteps
#                          
#      Functions to create:   Slope --> ouputs the slope values of density profile
#  

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
    path='/home/jcalahan/HighRes/'+filename

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
    
    # store distances of particles from COM (array)
    dist=np.sqrt((x-COMP[0])**2+(y-COMP[1])**2+(z-COMP[2])**2)
    
    #Initilize total mass
    mass=0
    
    # 1 particle of dark matter mass (assuming all dark matter particles have the same mass)
    dm_mass = m[0]
    
    #if a particle is within the limit, add a dark matter mass to total mass
    for d in dist:
        if d <= limit:
            mass = mass + dm_mass
                  
    return mass*(1e10)

def Density(snap,limit,steps):
    #Inputs: snap -- time step we want to study M33 at
    #        limit -- distance in kpc from the COM that defines the core of M33
    #        steps -- number of steps to measure density at
    #
    #Outputs: density -- array of floats, total dark matter mass/circular volume within limit at each step
    #         radii -- array of floats, radii at which we found the density
    #---------------------------------------
    
    #Define the step size 
    stepsize=float(limit)/float(steps)
    
    #Define the initial array of density values
    density=np.zeros(steps-1)
    radii=np.zeros(steps-1)
    
    #Initilize location and counter
    loc=stepsize
    n=0
    
    #Define Volume constants
    C=(4./3.)*np.pi
    
    #Find the density at intervals defined by stepsize
    while loc < limit:
        density[n]=ReadMass(snap,loc)/(C*loc**3)
        radii[n]=loc
        n+=1
        loc+=stepsize
        
    return density,radii

def func(r,a,c,d):
    return c*r**a

def Slope(r1,r2,density,radii):
    #Inputs: r1,r2 -- limits, find the slope between r1 and r2, r1 < r2
    #        density -- an array of density values as calculated from Density
    #        radii -- an array of radii values as calculated from Density (same size as 'density' array)
    #
    #Outputs: slope -- slope value from a linar fit
    #         
    #---------------------------------------
    
    #Initialize lists that will contain the sub group of radii and densities that we want to find the slope of
    sub_radii=[]
    sub_density=[]
    
    #Find radii and density that are within r1 and r2
    for i in range(len(density)):
        if radii[i] >= r1 and radii[i] <= r2:
            sub_radii.append(radii[i])
            sub_density.append(density[i])
            
    popt,pcov=curve_fit(func, sub_radii, sub_density)
    
    return popt
    
    
def ArrayofSlopes(s1,s2,numb,limit,steps): 
    #Inputs: s1,s2 -- two snap numbers to find the power-law fluxuations of between
    #        numb -- number of power-laws to calculate between s1 and s2
    #        limit -- distance in kpc from the COM that defines the core of M33
    #        steps -- number of steps to measure density at
    #
    #Outputs: slopes -- array of slopes
    #         time -- array of times that correspond to each slope
    #         
    #---------------------------------------
    
    #Initilize array of power-law values and time values
    slopes=np.zeros(numb+1)
    time=np.zeros(numb+1)
    snaps=[]
    
    #time intervals to examine the power-law value
    step=(s2-s1)/numb
    
    #Create a list of snap numbers
    s=s1
    while s<=s2:
        snaps.append(s)
        s=int(step+s)
    
    #Create array of times from snap numbers
    snapsarray=np.asarray(snaps)
    time=(snapsarray*3)/200
    #Initilize a counter
    i=0
    
    #Find the power-law values at snap numbers between s1 and s2
    while i <= numb:
        density,radii = Density(snaps[i],limit,steps)
        a=Slope(0.6,limit,density,radii)
        slopes[i]=a[0]
        #time[i]=(snaps[i]*3)/200
        i+=1
        
    return slopes,time
    



    
    