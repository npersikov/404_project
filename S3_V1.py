import numpy as np
import matplotlib.pyplot as plt

plt.close('all') # close plots

# Physical Parameters (those pertaining to droplets and forces)
dropMinD        = 0.0001; # Dropplet Min Diamater [m]
dropMaxD        = 0.1; # Dropplet Max Diamater [m]
initial_D       = np.random.random()*dropMaxD + dropMinD;
initial_h       = 0.1; # Initial Height [m]
rhoD            = 1000; # Density of dropplets [kg/m**3]
rhoA            = 1.225; # Air Density [kg/m**3]
cD              = .04; # Drag Coefficient for streamlined body [-]
surf_tension    = 0.25;
x0              = 0;

v0x         = 0.1; # Initial Velocity X [m/s]
v0y         = 0; # Initial Velocity Y [m/s]

def sim(d, h, x0, rhoD, rhoA, cD, v0x, v0y, surf_tension): # Look I turned the simulation into a function
    
    # Simulation Parameters (related to the numerical integration method)
    dt              = .001; # Time Step [s]
    # I               = 1; # Initial size of arrays for the simulation. This may be increased later.
    # size_increment  = 1; # how many indeces to increase array length by, if needed
    
    # Initialize droplet drag component arrays
    forceDragX = [];
    forceDragY = [];
    
    # Initialize X position and velocity arrays
    X = [x0];
    vX = [v0x];
    
    # Initialize Y position and velocity arrays
    Y = [h];
    vY = [v0y];
    
    area = 0.25*np.pi*d**2; # Calculate areas
    mass = rhoD*(4/3)*np.pi*(d/2)**3; # Generate masses
    
    # Main Simulation Loop
    i = 0; # start the position index at 0 for each droplet at the beginning
    
    weber_number = 0;
    while weber_number <= 1: # propagate until it breaks up. TODO change this condition!!!
        
        # if the loop must continue past the initially allocated size, resize the arrays and update their length.
        # if i >= array_size - 1:   
        #     X.resize(len(X) + size_increment);
        #     Y.resize(len(Y) + size_increment);
        #     vX.resize(len(vX) + size_increment);
        #     vY.resize(len(vY) + size_increment);
        #     forceDragX.resize(len(forceDragX) + size_increment);
        #     forceDragY.resize(len(forceDragY) + size_increment);
        #     array_size += size_increment; 
    
        # Calculate drag force components
        forceDragX = np.append(forceDragX, -np.sign(vX[i])*0.5*rhoA*cD*area*vX[i]**2); 
        forceDragY = np.append(forceDragY, -np.sign(vY[i])*0.5*rhoA*cD*area*vY[i]**2); #should have different sign from x
        
        # Sum the forces and solve for acceleration
        aX = forceDragX[i]/mass;
        aY = forceDragY[i]/mass - 9.81;
        
        # Integrate acceleration to get velocity
        vX = np.append(vX, vX[i]+aX*dt);
        vY = np.append(vY, vY[i]+aY*dt);
    
        # Integrate velocity to get position
        X = np.append(X, X[i]+vX[i]*dt);
        Y = np.append(Y, Y[i]+vY[i]*dt);
            
        i += 1; # Increment the simulation index
        weber_number = rhoD*((vX[i]**2)+(vY[i]**2))*d/surf_tension;
    
    return X, Y, vX, vY; # End of the simulation function

def recursive_mission_run(d, h, x0, rhoD, rhoA, cD, v0x, v0y, surf_tension):
    X, Y, vX, vY = sim(d, h, x0, rhoD, rhoA, cD, v0x, v0y, surf_tension);
    plt.plot(X,Y);
    if Y[len(Y) - 1] > 0:
        # print("recursive run")
        recursive_mission_run(d/2, Y[len(Y) - 1], X[len(X) - 1], rhoD, rhoA, cD, vX[len(vX) - 1]*1.1, vY[len(vX) - 1]*0.9, surf_tension)
        recursive_mission_run(d/2, Y[len(Y) - 1], X[len(X) - 1], rhoD, rhoA, cD, vX[len(vX) - 1]*0.9, vY[len(vX) - 1]*1.1, surf_tension)
    else:
        return;   

h = initial_h;
d = initial_D;
paths = plt.figure()
plt.ylabel('Y Position (m)')
plt.xlabel('X Position (m)')
plt.title('Droplet Flight Paths')
recursive_mission_run(d, h, x0, rhoD, rhoA, cD, v0x, v0y, surf_tension);



# while h > 0:
#     X, Y, vX, vY = sim(d, h, x0, rhoD, rhoA, cD, v0x, v0y, surf_tension);
#     zero_heights, = np.where(Y[0:len(Y)] <= 0);
#     ground_hit_index = zero_heights[0];
#     plt.plot(X[0:ground_hit_index],Y[0:ground_hit_index]);
#     d = d/2;
#     h = Y[len(Y) - 1];
#     v0x = 
    

