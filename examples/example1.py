#!/usr/bin/python

import pyvicsek as vi

# We define the simulation parameters

N_steps = 5 # Number of steps per measure
N_it = 400 # Number of measures

N_birds = 300 # Number of birds
length = 7 # Length of the box
vel = .03 # Speed of the birds
ang_noise = 2. # Angular noise of the copying mechanism
SEED = 2302 # Seed for the random number generation


# We create the simulation object
sim = vi.VicsekSimulation(N_birds, L=length, v=vel, eta=ang_noise, seed=SEED)

with open("example_va.txt", "w") as f:
    f.write("#t\tv_a\n") # We write a header for the file
    f.write(f"{sim.t}\t{sim.compute_va()}\n") # We write the initial conditions

    for _ in range(N_it): # Simulation loop
        sim.update(N_steps) # We update the simulation
        f.write(f"{sim.t}\t{sim.compute_va()}\n") # We write the current order parameter
