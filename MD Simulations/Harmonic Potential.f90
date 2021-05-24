### alternate initialization of constants are usedin case of sho



import math
import matplotlib.pyplot as plt
import numpy as np
mass =  145123   ## or 1837
omega = 0.01
sigma = 6.69
dt = 75             #1/50 * 2* (math.pi) /omega
total_time = 10**6  #20*2*(math.pi)/omega
number_steps = int(total_time/dt)
x=7
v=0
#x0 = 7.5
eps = 8.71*(10**(-4))
current_time = 0


def evolve_by_1step(x,v,acc):
    x += v*dt + 0.5*acc*dt*dt               
    v += 0.5*acc*dt                         # using velocity- verlet
    acc,pot = calculate_acceleration(x)
    v += 0.5*acc*dt
    energy = pot + 0.5*mass*(v**2)
    return x,v,acc,pot,energy
    
def calculate_acceleration(x):
    #acc= -1*(omega**2)*(x - x0)            
    acc = 4*eps*(12*((sigma**12)/ x**(13))-6*((sigma**6)/(x**7))) / mass   #(a = -1/m *(dV/dx)) 
    #pot = 0.5*mass*(omega**2)*(x-x0)**2    
    pot = 4*eps*((sigma/x)**12 - (sigma/x)**6)
    return acc,pot

acc,pot = calculate_acceleration(x)

file_coords = open("coordinates.out","w")
energy_file = open("energy.out","w")
potfile     = open("potf.out","w")
current_time = 0

for i in range(number_steps):  
    print(current_time,x,v,pot,file=file_coords)            #(printing coords vs time)
    x,v,acc,pot,energy = evolve_by_1step(x,v,acc)
    print(current_time,energy,file=energy_file)
    current_time = current_time + dt
    

for i in range(int(number_steps)):  
    x,v,acc,pot,energy = evolve_by_1step(x,v,acc)
    print(x,pot,file=potfile)
    current_time = current_time + dt



file_coords.close()
energy_file.close()
potfile.close()

data = np.loadtxt("coordinates.out")
plt.plot(data[:,0],data[:,1]) 
plt.xlabel("Time") # Text for X-Axis
plt.ylabel("Position") # Text for Y-Axis
plt.show()   

data = np.loadtxt("potf.out")
plt.plot(data[:,0],data[:,1]) 
plt.show()  

data_energy = np.loadtxt("energy.out")
plt.plot(data_energy[:,0],data_energy[:,1])   
plt.show() 
