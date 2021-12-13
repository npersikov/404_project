
import numpy as np
import matplotlib.pyplot as plt

plt.close('all') # close plots

# Physical Parameters (those pertaining to droplets and forces)
dropMinD    = 0.0001; # Dropplet Min Diamater [m]
dropMaxD    = 0.1; # Dropplet Max Diamater [m]
numDrops    = 10; # Number of Dropplets [n]
h           = 500; # Initial Height [m]
rhoD        = 1000; # Density of dropplets [kg/m**3]
v0x         = 1000; # Initial Velocity X [m/s]
v0y         = 0; # Initial Velocity Y [m/s]
rhoA        = 1.225; # Air Density [kg/m**3]
cD          = .04; # Drag Coefficient for streamlined body [-]

# Simulation Parameters (related to the numerical integration method)
dt              = .001; # Time Step [s]
I               = 1000; # Initial size of arrays for the simulation. This may be increased later.
size_increment  = 100; # how many indeces to increase array length by, if needed

# Initialize droplet diameter and area arrays
dropDs = np.zeros(numDrops);
dropAs = np.zeros(numDrops);

# Initialize droplet mass and drag component arrays
massD = np.zeros(numDrops);
forceDragX = np.zeros((I, numDrops));
forceDragY = np.zeros((I, numDrops));

# Initialize X position and velocity arrays
X = np.zeros((I, numDrops));
vX = np.zeros((I, numDrops));
vX[0,:] = v0x;

# Initialize Y position and velocity arrays
Y = np.zeros((I, numDrops));
Y[0,:] = h; # set first y positions of each particle to h
vY = np.zeros((I, numDrops));

# Generate droplets with diameters and masses
for i in range(numDrops) :
    dropDs[i] = np.random.random()*dropMaxD + dropMinD; # Generate diameters
    dropAs[i] = 0.25*np.pi*dropDs[i]**2; # Calculate areas
    massD[i] = rhoD*(4/3)*np.pi*(dropDs[i]/2)**3; # Generate masses

# Main Simulation Loop
for j in range(numDrops) : # For each droplet
    i = 0; # start the position index at 0 for each droplet at the beginning
    array_size = I; # set array size to the initial allocation, which is 1000
    while Y[i,j] > 0: # propagate until it hits the ground
        
        # if the loop must continue past the initially allocated size, resize the arrays and update their length.
        if i >= array_size - 1:   
            X.resize(len(X) + size_increment, numDrops);
            Y.resize(len(Y) + size_increment, numDrops);
            vX.resize(len(vX) + size_increment, numDrops);
            vY.resize(len(vY) + size_increment, numDrops);
            forceDragX.resize(len(forceDragX) + size_increment, numDrops);
            forceDragY.resize(len(forceDragY) + size_increment, numDrops);
            array_size += size_increment; 
    
        # Calculate drag force components
        forceDragX[i,j] = -np.sign(vX[i,j])*0.5*rhoA*cD*dropAs[j]*vX[i,j]**2; 
        forceDragY[i,j] = -np.sign(vY[i,j])*0.5*rhoA*cD*dropAs[j]*vY[i,j]**2; #should have different sign from x
        
        # Sum the forces and solve for acceleration
        aX = forceDragX[i,j]/massD[j];
        aY = forceDragY[i,j]/massD[j] - 9.81;
        
        if abs(aY) < 0.01 and i % 100 == 0:
            print("Droplet ", j , " has reached terminal velocity!");
        
        if i < array_size - 1:
            # Integrate acceleration to get velocity
            vX[i+1,j] = vX[i,j]+aX*dt;
            vY[i+1,j] = vY[i,j]+aY*dt;
        
            # Integrate velocity to get position
            X[i+1,j] = X[i,j]+vX[i,j]*dt;
            Y[i+1,j] = Y[i,j]+vY[i,j]*dt;
            
        i += 1; # Increment the simulation index
        
# Plot the droplet paths. Each drop must be plotted separately since they 
# travel for a different amount of time. The index at which the drop hits the 
# ground is found using the "np.where" function.
paths = plt.figure()
plt.ylabel('Y Position (m)')
plt.xlabel('X Position (m)')
plt.title('Droplet Flight Paths')
for drop in range(numDrops):
    zero_heights, = np.where(Y[0:len(Y), drop] <= 0);
    ground_hit_index = zero_heights[0];
    plt.plot(X[0:ground_hit_index, drop],Y[0:ground_hit_index, drop]);

# Plot the masses of the droplets against their impact point
masses = plt.figure()
plt.ylabel('Droplet Landing Point (m)')
plt.xlabel('Droplet Mass (kg)')
plt.title('Droplet Flight Distance vs Droplet Mass')
for drop in range(numDrops):
    zero_heights, = np.where(Y[0:len(Y), drop] <= 0);
    ground_hit_index = zero_heights[0];
    plt.scatter(massD[drop], X[ground_hit_index, drop]);
