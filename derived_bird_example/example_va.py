#!/usr/bin/python

import pyvicsekvision as viv

import numpy as np # np.pi, np.random.rand

# We define the simulation parameters

N_steps = 5 # Number of steps per measure
N_it = 400 # Number of measures

N_birds = 300 # Number of birds
length = 7 # Length of the box
vel = .03 # Speed of the birds
ang_noise = 2. # Angular noise of the copying mechanism
width_of_view = np.pi/3 # Wdith of the field of view of the birds
SEED = 2302 # Seed for the random number generation

# We create a bird list. In order to use our new kind of birds a list is required
birdlist = [viv.BirdVision(width_of_view) for _ in range(N_birds)]

# We create the simulation object
sim = viv.VicsekSimulation(birdlist, L=length, v=vel, eta=ang_noise, seed=SEED)

# We set the coordinates of the birds to random, this is required
# since in the list declaration we set them to the default value
# which is (x, y, theta) = (0, 0, 0)
sim.randomize()

with open("example_va.txt", "w") as f:
    f.write("#t\tv_a\n") # We write a header for the file
    f.write(f"{sim.t}\t{sim.compute_va()}\n") # We write the initial conditions

    for _ in range(N_it): # Simulation loop
        sim.update(N_steps) # We update the simulation
        f.write(f"{sim.t}\t{sim.compute_va()}\n") # We write the current order parameter
