# This code simulates the attomization of 
# fluid exiting a pipe of fixed diameter, headed 
# at initial velocity in a downward trajectory. 
# As drag force overcomes cohesion force
# the fluid breaks up, creating smaller and smaller
# dropplets until drag force no longer can excede 
# cohesion and the dropplets reach terminal velocity.


# Expected output plots:
# Time v Dropplet size (to show the breakup)
# velocity v Y position (to show where terminal velocity is reached)(also shows breakup radius)

import numpy as np
import matplotlib.pyplot as plt
######################################################### KEY VARIABLES
# Pipe Diameter [m]
dP = 0.01

# Initial Velocity Y [m/s] 
# Note: negative because pipe pointed downward
v0y = -1

# Air Density [kg/m**3]
rhoA = 1.225

# Drag Coefficient for streamlined body [-]
cD = .04

# Fluid density [kg/m**3]
rhoF = 1000

# Time Step [s]
dt = .001

# Iterations [n]
I = 50000

#acelleration due to gravity [m/s**2]
g = -9.81

# how low the acceleration has to get to declare 
# the dropplet has reached terminal velocity [m/s**2]
margin = -.00001


################################################### ARRAY INITIALIZATION

# initialize arrays. will be appended to until calculations converge
# End condition (acceleration ~= 0)

Y = [0]
vY = [v0y]

# drag and acceleration both do not have initial values since 
# they require one loop of calculation to compute the first value
forceDragY = []
aY = []

# time array
T = [0]




############################################### DEPENDANT VARIABLES
currentAcc = g

# drop diameter is assumed to be the same as the pipe
dropD = dP

# Generate areas
dropA = 0.25*np.pi*dropD**2

# Generate masses of dropplets
massD = rhoF*(4/3)*np.pi*(dropD/2)**3


############################################### ITTERATIVE CALCS
# Position Generator
t = 0
i = 0
while currentAcc < margin:
    
   

    #drag force in Y
    forceDragY = np.append(forceDragY, 0.5*rhoA*cD*dropA*vY[i]**2) #should have different sign from x
    
    # NOTE drag is assumed to be negative for x and positive for y. This may not be the case. Ideally, it would depend on the velocity vector
    # aY = 0.5*rhoA*cD*dropAs[j]*vY[i,j]**2/massD[j] - 9.81;
    aY = np.append(aY, forceDragY[i]/massD + g)
    currentAcc = aY[i]
      
    # update velocity
    vY = np.append(vY,vY[i]+aY[i]*dt)
    
    # Update position
    Y = np.append(Y, Y[i]+vY[i]*dt)
    
    # increase iteration step by 1
    i = i+1 
    
    # Increase time step
    t = t+dt
    T = np.append(T,t)
    
############################################## PLOTS    
        # Plot for dropplets
        
# Position
paths = plt.figure()
#plt.ylim(0,110)
plt.ylabel('Y Position (m)')
plt.xlabel('Time (s)')
plt.title('Droplet Position')
plt.plot(T,Y)

# Velocity
paths = plt.figure()
#plt.ylim(0,110)
plt.ylabel('Y Velocity (m/s)')
plt.xlabel('Time (s)')
plt.title('Droplet Velocity')
plt.plot(T,vY)

# Acceleration
paths = plt.figure()
#plt.ylim(0,110)
plt.ylabel('Y Acceleration (m/s^2)')
plt.xlabel('Time (s)')
plt.title('Droplet Acceleration')
plt.plot(T[1:],aY)

''' saving structure for plot of dropplet size v time
   #Plot for mass v distance
masses = plt.figure()
plt.ylabel('Droplet Landing Point (m)')
plt.xlabel('Droplet Mass (kg)')
plt.title('Droplet Flight Distance vs Droplet Mass')
for j in range(numDrops) :
    plt.scatter(massD, X[len(X)-1])
    '''
    
    
 ############################################### PRINT RESULTS   
print("time to terminal velocity: " + str(max(T)) + " seconds" )



#####S2_V1
