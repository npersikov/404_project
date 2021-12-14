import numpy as np
import matplotlib.pyplot as plt

plt.close('all') # close plots

initial_D       = 0.01;
initial_h       = 10; # Initial Height [m]
rhoD            = 1000; # Density of dropplets [kg/m**3]
rhoA            = 1.225; # Air Density [kg/m**3]
cD              = .04; # Drag Coefficient for streamlined body [-]
surf_tension    = 0.072;
x0              = 0;

v0x         = 0.2; # Initial Velocity X [m/s]
v0y         = 0; # Initial Velocity Y [m/s]

def sim(d, h, x0, rhoD, rhoA, cD, v0x, v0y, surf_tension): # Look I turned the simulation into a function
    
    # Simulation Parameters (related to the numerical integration method)
    dt              = .01; # Time Step [s]
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
    while weber_number <= 11 and Y[i] > 0: # propagate until it breaks up. TODO change this condition!!!
    
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
        weber_number = rhoA*((vX[i]**2)+(vY[i]**2))*d/surf_tension;
    
    return X, Y, vX, vY; # End of the simulation function

def get_mass(d, rhoD):
    volume = 4/3*np.pi*(d/2)**3
    return rhoD*volume

def get_diameter(rhoD, mass):
    volume = mass/rhoD;
    d = 2*(3/4*volume/np.pi)**(1/3);
    return d;

def recursive_mission_run(d, h, x0, rhoD, rhoA, cD, v0x, v0y, surf_tension, paths, masses):
    X, Y, vX, vY = sim(d, h, x0, rhoD, rhoA, cD, v0x, v0y, surf_tension);
    paths.plot(X,Y);
    xmin, xmax = paths.get_xlim();
    if xmax < X[len(X) - 1]:
        paths.set_xlim([0, X[len(X) - 1]*1.1]);
    if Y[len(Y) - 1] > 0:
        # print("recursive run")
        
        min_v_reduction_factor = 0.85;
        max_v_reduction_factor = 0.95;
        
        new_mass_1  = get_mass(d, rhoD)*(0.25+0.5*np.random.random())
        new_mass_2  = get_mass(d, rhoD)-new_mass_1
        new_xV_1    = vX[len(vX) - 1]*(min_v_reduction_factor + (max_v_reduction_factor - min_v_reduction_factor)*np.random.random())
        new_yV_1    = vY[len(vY) - 1]*(min_v_reduction_factor + (max_v_reduction_factor - min_v_reduction_factor)*np.random.random())
        new_x_mom_1 = new_mass_1*new_xV_1
        new_y_mom_1 = new_mass_1*new_yV_1
        
        old_momentum_x = get_mass(d, rhoD) * vX[len(vX) - 1];
        old_momentum_y = get_mass(d, rhoD) * vY[len(vY) - 1];
        
        new_x_mom_2 = old_momentum_x - new_x_mom_1
        new_y_mom_2 = old_momentum_y - new_y_mom_1
        new_xV_2    = new_x_mom_2/new_mass_2
        new_yV_2    = new_y_mom_2/new_mass_2
        # new_d_1     = 2*(3*new_mass_1*rhoD/(4*np.pi))**(1/3)
        # new_d_2     = 2*(3*new_mass_2*rhoD/(4*np.pi))**(1/3)
        new_d_1     = get_diameter(rhoD, new_mass_1)
        new_d_2     = get_diameter(rhoD, new_mass_2)
        
        recursive_mission_run(new_d_1, Y[len(Y) - 1], X[len(X) - 1], rhoD, rhoA, cD, new_xV_1, new_yV_1, surf_tension, paths, masses)
        recursive_mission_run(new_d_2, Y[len(Y) - 1], X[len(X) - 1], rhoD, rhoA, cD, new_xV_2, new_yV_2, surf_tension, paths, masses)
    else:
        zero_heights, = np.where(Y[0:len(Y)] <= 0);
        ground_hit_index = zero_heights[0];
        
        x_before = X[ground_hit_index - 1]
        x_after = X[ground_hit_index]
        y_before = X[ground_hit_index - 1]
        y_after = X[ground_hit_index]
        
        # plot the mass at the point closest to y = 0 index
        if abs(y_before) > abs(y_after):
            # masses.scatter(get_mass(d, rhoD), x_after);
            masses.scatter(x_after, get_mass(d, rhoD));
        else:
            # masses.scatter(get_mass(d, rhoD), x_before);
            masses.scatter(x_before, get_mass(d, rhoD));
        
        return;   

h = initial_h;
d = initial_D;
# paths = plt.figure()
# paths.ylabel('Y Position (m)')
# paths.xlabel('X Position (m)')
# paths.title('Droplet Flight Paths')
# paths.ylim([0, initial_h*1.1])
# paths.xlim([0, 0.05]);
paths = plt.subplot(1,2,1)
paths.set_ylabel('Y Position (m)')
paths.set_xlabel('X Position (m)')
paths.set_title('Droplet Flight Paths')
paths.set_ylim([0, initial_h*1.1])
paths.set_xlim([0, v0x*np.sqrt(2*initial_h/9.81)*1.1]);

masses = plt.subplot(1,2,2)
masses.set_ylabel('Mass (kg)')
masses.set_xlabel('X Position (m)')
masses.set_title('Distance Travelled by Droplet Mass')
# masses.set_ylim([0, initial_h*1.1])
# masses.set_xlim([0, 0.05]);

recursive_mission_run(d, h, x0, rhoD, rhoA, cD, v0x, v0y, surf_tension, paths, masses);
    

