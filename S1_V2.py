
import numpy as np
import matplotlib.pyplot as plt

dropMinD = 0.0001;

# Dropplet Max Diamater [m]
dropMaxD = 0.1;

# Number of Dropplets [n]
numDrops = 10;

# Initial Height [m]
h = 100;

# Density of dropplets [kg/m**3]
rhoD = 1000;

# Initial Velocity X [m/s]
v0x = 1000;

# Initial Velocity Y [m/s]
v0y = 0;

# Air Density [kg/m**3]
rhoA = 1.225;

# Drag Coefficient for streamlined body [-]
cD = .04;

# Time Step [s]
dt = .001;

# Initial size of arrays for the simulation. This may be increased later.
I = 1000;

dropDs = np.zeros(numDrops);
dropAs = np.zeros(numDrops);

massD = np.zeros(numDrops);
forceDragX = np.zeros((I, numDrops));
forceDragY = np.zeros((I, numDrops));

X = np.zeros((I, numDrops));
vX = np.zeros((I, numDrops));
vX[0,:] = v0x;
# aX = np.zeros((I, numDrops));

Y = np.zeros((I, numDrops));
Y[0,:] = h; # set first y positions of each particle to h
vY = np.zeros((I, numDrops));
# aY = np.zeros((I, numDrops));

# Calculations
for i in range(numDrops) :
    # Generate Dropplet Diameters within range
    D = 0;
    
    while D < dropMinD :
        D = np.random.random()*dropMaxD;
    
    dropDs[i] = D;

    # Generate areas
    dropAs[i] = 0.25*np.pi*dropDs[i]**2;

    # Generate masses of dropplets
    massD[i] = rhoD*(4/3)*np.pi*(dropDs[i]/2)**3;

# Position Generator
t = 0;
i = 0;
size_increment = 100; # how many indeces to increase array length by, if needed

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
            # aX.resize(len(aX) + size_increment, numDrops);
            # aY.resize(len(aY) + size_increment, numDrops);
            forceDragX.resize(len(forceDragX) + size_increment, numDrops);
            forceDragY.resize(len(forceDragY) + size_increment, numDrops);
            array_size += size_increment; 
            #print("resizing arrays to: ", array_size); 
    
        # calculate drag forces
        forceDragX[i,j] = -0.5*rhoA*cD*dropAs[j]*vX[i,j]**2; 
        forceDragY[i,j] = 0.5*rhoA*cD*dropAs[j]*vY[i,j]**2; #should have different sign from x
        
        # NOTE drag is assumed to be negative for x and positive for y. This may not be the case. Ideally, it would depend on the velocity vector
        aX = forceDragX[i,j]/massD[j];
        # aY = 0.5*rhoA*cD*dropAs[j]*vY[i,j]**2/massD[j] - 9.81;
        aY = forceDragY[i,j]/massD[j] - 9.81;
        
        if i < array_size - 1:
            # Apply drag
            
            vX[i+1,j] = vX[i,j]+aX*dt;
            vY[i+1,j] = vY[i,j]+aY*dt;
        
            # Generate X and Y positions with drag and gravity
            X[i+1,j] = X[i,j]+vX[i,j]*dt;
            Y[i+1,j] = Y[i,j]+vY[i,j]*dt;
            
        i += 1;            
        
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
    
paths2 = plt.figure()
plt.ylabel('Y Position (m)')
plt.xlabel('X Position (m)')
plt.title('Droplet Flight Paths 2')
plt.plot(X,Y);
# for drop in range(numDrops):
#     zero_heights, = np.where(Y[0:len(Y), drop] <= 0);
#     ground_hit_index = zero_heights[0];
#     plt.plot(X[0:ground_hit_index, drop],Y[0:ground_hit_index, drop]);

masses = plt.figure()
plt.ylabel('Droplet Landing Point (m)')
plt.xlabel('Droplet Mass (kg)')
plt.title('Droplet Flight Distance vs Droplet Mass')
for drop in range(numDrops):
    zero_heights, = np.where(Y[0:len(Y), drop] <= 0);
    ground_hit_index = zero_heights[0];
    plt.scatter(massD[drop], X[ground_hit_index, drop]);
