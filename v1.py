
import numpy as np
import matplotlib.pyplot as plt

dropMinD = 0.0001;

# Dropplet Max Diamater [m]
dropMaxD = 0.1;

# Number of Dropplets [n]
numDrops = 10;

# Initial Height [m]
h = 1;

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

# Iterations [n]
I = 10000;

dropDs = np.zeros(numDrops);
dropAs = np.zeros(numDrops);

massD = np.zeros(numDrops);
forceDrag = np.zeros((I, numDrops));

X = np.zeros((I, numDrops));
vX = np.zeros((I, numDrops));
vX[0,:] = v0x;
aX = np.zeros((I, numDrops));

Y = np.zeros((I, numDrops));
Y[0,:] = h; # set first y positions of each particle to h
vY = np.zeros((I, numDrops));
aY = np.zeros((I, numDrops));

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
for i in range(I) :
    t = t+dt;
    for j in range(numDrops) :
        
        # TODO do something with this, rename for x and make another for y, get rid of it, or save the magnitude
        forceDrag[i,j] = -0.5*rhoA*cD*dropAs[j]*vX[i,j]**2; 
        
        # NOTE drag is assumed to be negative for x and positive for y. This may not be the case. Ideally, it would depend on the velocity vector
        aX = forceDrag[i,j]/massD[j];
        # aY = 0.5*rhoA*cD*dropAs[j]*vY[i,j]**2/massD[j] - 9.81;
        aY = 0.5*rhoA*cD*dropAs[j]*vY[i,j]**2/massD[j] - 9.81;
        
        if i < I-1:
            # Apply drag
            
            vX[i+1,j] = vX[i,j]+aX*dt;
            vY[i+1,j] = vY[i,j]+aY*dt;
        
            # Generate X and Y positions with drag and gravity
            X[i+1,j] = X[i,j]+vX[i,j];
            Y[i+1,j] = Y[i,j]+vY[i,j];
        

for i in range(numDrops) :
    plt.plot(X,Y);