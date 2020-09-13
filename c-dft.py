import numpy as np

#ideal noninteracting fluid.
def ideal_p(distance, density):
    i_fluid = 1/beta * distance * density * (np.log((delta**3)*density) - 1)
    return i_fluid
#Lenard-Jones potential
def lenard_jones(d):
    #Epsilon = 1.6 KJ/mol
    #sigma = 2.6 Angstroms
    lj_potential = 4 * 1.6 * ( ( (2.6/d)**12 ) - ((2.6/d)**6) )
    return lj_potential

# dummy_correlation
def dummy_correlation(prev_density,calc_density):
    diff_density = calc_density - prev_density
    if diff_density < 0:
        new_density =  (calc_density + (diff_density)) * -0.05
    else:
        new_density =  (calc_density + (diff_density)) * 0.05
    return new_density

#Variables
distances = np.linspace(1,10,100)# angstroms, starting distance, ending distance, points between the start and the end
densities = np.linspace(0.01,0.1,100)# molecules/angstrom^3 starting density, ending density, points between the start and the end

# Constants
kb = 1.38065E-26 # Boltzmann constant, KJ/K
T = 300 # Temparute, K 
Na = 6.022E23 # Avogrado's constant, 1/mol
atom_mass = 20 # mass, u
plank_constant = 6.26E-37 # KJ * s
bulk_density = 0.033 # molecules / angstrom^3 

#Calculated constants
beta = 1 / (Na*kb*T) # Boltzmann factor, KJ / mol
delta = np.sqrt(beta*(plank_constant**2)/(2*np.pi*atom_mass)) #
chem_potential = np.log(delta**3 * bulk_density)/beta # 

#is used in the analytical solution
a_distance_v_density = {}
#is used in the numerical aproximation

distance_v_density = {} # will hold all the disctances as keys and the densities (the ones that lowered the total energy)
prev_density = 0 
for i in range(len(distances)-1):
    helmoltz_v_density = {} #will hold the total energy per density as key and the calculated density as value
    distance = distances[i] # get a distance from the list of distances
    ext_p = lenard_jones(distance) # define an external potential here
    
    #analytical solution
    calc_density = bulk_density * np.exp(-beta*ext_p)
    correlation = dummy_correlation(prev_density,calc_density)
    a_distance_v_density[distance] = calc_density + correlation
    prev_density = calc_density + correlation
    #Numerical aproximation
    # calculate the total energy for every density in the list
    for density in densities:
         
        e_tot = ideal_p(distance,density) + (distance*ext_p*density) - distance*chem_potential*density
        helmoltz_v_density[e_tot] = density
    #search for the minimum total energy and added as a value of the given distance
    distance_v_density[distance] = helmoltz_v_density[min(helmoltz_v_density)]
    print(f"{distance},{distance_v_density[distance]},{a_distance_v_density[distance]}") #simple output



