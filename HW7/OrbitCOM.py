# Homework 6 Solution
# G. Besla 


# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# NOTE: I have modified CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX
from CenterOfMass import CenterOfMass




def OrbitCOM(galaxy,start,end,n):
    # function that loops over all the desired snapshots to compute the COM pos and vel
    # as a function of time
    # inputs
    #       galaxy  the name of the galaxy e.g. "MW"
    #       start   initial SnapNumber  e.g. 0
    #       end     final SnapNumber  e.g. 100
    #       n       integer indicating the frequency with which the snapshots are read in;  n should not be 0
    # creates:  file with t, x, y, z, vx, vy, vz for n snapshots
      
    
    fileout = 'Orbit_%s.txt'%(galaxy)   # filename to store output
      
    # set tolerance for COM convergence
    delta = 1.0
    # decide by how much will decrease RMAX
    VolDec = 2.0
    # For M33 there needs to be a smaller tolerance, as the object gets stripped the centroiding gets worse
    if (galaxy == "M33"):
        delta = 0.2
        VolDec = 4.0
    
    a = int(end/n)+1
    Orbit = np.zeros((a,7))         # initialize array that will store the COM position
        # np.zeros initializes with 0s 
        # in the end we want: t, xC, yC, zC, vxC, vyC, vzC  = 7 columns
        # row, column
        
    for i in np.arange(start, end+n, n):            # loop over the files in intervals of n

        # Determine Filename
        ilbl = '000' + str(i)                       # add a string of the filenumber to the value "000"
        ilbl = ilbl[-3:]                            # remove all but the last 3 digits
        filename='%s_'%(galaxy) + ilbl + '.txt'     # create filenames
            
                      
        COM = CenterOfMass(filename,2)              # Create a Center of Mass Object using Disk Particles
          
        Orbit[int(i/n),0] = float(COM.time/u.Myr/1000.)    # Store time in Gyr                  
            # time is stored in the CenterOfMass object 
            # when you called Read in the CenterOfMass Class you initialized a global variable self.time
            # when you create a new instance of this class, self.time is now a global variable 
            # associated with the object
            # need to get rid of the units when storing in array - conver the scalar "quantity" to a float
            
        COMP = COM.COM_P(delta,VolDec)   # Define the COM position
        COMV = COM.COM_V(COMP[0],COMP[1],COMP[2])  # Define the COM velocity
            
        #Need to get rid of units when storing in array and then convert the scalar "quantity" to a float
        Orbit[int(i/n),1] = float(COMP[0]/u.kpc)  # Store COM Position. 
        Orbit[int(i/n),2] = float(COMP[1]/u.kpc)
        Orbit[int(i/n),3] = float(COMP[2]/u.kpc)
        Orbit[int(i/n),4] = float(COMV[0]/u.km*u.s) # Store COM Velocity
        Orbit[int(i/n),5] = float(COMV[1]/u.km*u.s)
        Orbit[int(i/n),6] = float(COMV[2]/u.km*u.s)
    

        print(i)  # print the counter so you know how far the code has progressed.
        
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, Orbit, header='t   x    y    z   vx   vy   vz', comments='# ',
                fmt=['%.2f', '%.2f','%.2f','%.2f','%.2f','%.2f','%.2f']) 
      

        



# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
print("Starting to Compute MW Orbit")
OrbitCOM("MW",0,800,5)
print("Starting to Compute M31 Orbit")
OrbitCOM("M31",0,800,5)
print("Starting to Compute M33 Orbit")
OrbitCOM("M33",0,800,5)



# Read in the data files
# headers:  t, x, y, z, vx, vy, vz

MWOrbit = np.genfromtxt('Orbit_MW.txt',dtype=None,names=True) 
M31Orbit = np.genfromtxt('Orbit_M31.txt',dtype=None,names=True) 
M33Orbit = np.genfromtxt('Orbit_M33.txt',dtype=None,names=True) 



# Determine the magnitude of the 
# relative position and velocities 

# of MW and M31
M31MWR = np.sqrt((MWOrbit['x']-M31Orbit['x'])**2 + (MWOrbit['y']-M31Orbit['y'])**2 
                 +(MWOrbit['z']-M31Orbit['z'])**2) 
M31MWV = np.sqrt((MWOrbit['vx']-M31Orbit['vx'])**2 + (MWOrbit['vy']-M31Orbit['vy'])**2 
                 +(MWOrbit['vz']-M31Orbit['vz'])**2) 

# of M33 and M31
M33M31R = np.sqrt((M31Orbit['x']-M33Orbit['x'])**2 + (M31Orbit['y']-M33Orbit['y'])**2 
                 +(M31Orbit['z']-M33Orbit['z'])**2) 
M33M31V = np.sqrt((M31Orbit['vx']-M33Orbit['vx'])**2 + (M31Orbit['vy']-M33Orbit['vy'])**2 
                 +(M31Orbit['vz']-M33Orbit['vz'])**2)




# Plot the Orbit of the galaxies 
#################################


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot the separtion of M31 and MW
plt.plot(MWOrbit['t'], M31MWR, color='blue', linewidth=5, label='MW-M31')

# Plot the separtion of M33 and M31
plt.plot(MWOrbit['t'], M33M31R, color='red', linewidth=5, label='M31-M33')

# Add axis labels
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Separation (kpc)', fontsize=22)

#set axis limits
plt.ylim(0,200)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')

# Save to a file
ax.set_rasterized(True)
plt.savefig('OrbitR.eps', rasterized=True, dpi=350)




# Plot a zoom of the MW-M31 Orbit 
#################################


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot the separtion of M31 and MW
plt.plot(MWOrbit['t'], M31MWR, color='blue', linewidth=5, label='MW-M31')

# Add axis labels
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Separation (kpc)', fontsize=22)

#set axis limits
plt.ylim(0,50)
plt.xlim(6,6.5)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')


# Save to a file
ax.set_rasterized(True)
plt.savefig('OrbitR_ZoomMWM31.eps', rasterized=True, dpi=350)




# Plot a zoom of the M33-M31 Orbit 
#################################


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot the separtion of M33 and M31
plt.plot(M33Orbit['t'], M33M31R, color='red', linewidth=5, label='M31-M33')


# Add axis labels
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Separation (kpc)', fontsize=22)

#set axis limits
plt.ylim(0,125)
plt.xlim(7,12)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')


# Save to a file
ax.set_rasterized(True)
plt.savefig('OrbitR_ZoomM33M31.eps', rasterized=True, dpi=350)



# Bonus Solution from Andy Henrici

# Estimation of M33 orbit decay after 6.0 Gyr
time = M33Orbit['t']
# Grab only the data after 6.0 Gyr
index = np.where(time > 6.0)
x = time[index]
y = M33M31R[index]

# define the exponential function to test fit
def model_func(x, a, k, b):
    return a * np.exp(-k*x) + b

# Plot the fit to the data
fig = plt.figure(2)
t = np.linspace(6.0, 20)
plt.plot(x,y,'o')
plt.xlabel('Time [Gyr]')
plt.ylabel('Distance [kpc]')
f = model_func(t, 300, 0.13, 0)
plt.plot(t, f, 'r-')




# Plot the Relative velocites of the galaxies 
####################################


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot the relative velocity of M31 and MW
plt.plot(MWOrbit['t'], M31MWV, color='blue', linewidth=5, label='MW-M31')

# Plot the relative velocity of M33 and M31
plt.plot(MWOrbit['t'], M33M31V, color='red', linewidth=5, label='M31-M33')

# Add axis labels
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Relative Velocity (km/s)', fontsize=22)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')


# Save to a file
ax.set_rasterized(True)
plt.savefig('OrbitV.eps', rasterized=True, dpi=350)





