#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import pyvicsekvision as viv

# We define the simulation parameters

N_steps = 5 # Number of steps between frames
N_frames = 400 # Number of frames

N_birds = 300 # Number of birds
length = 7 # Length of the box
vel = .03 # Speed of the birds
ang_noise = 2. # Angular noise of the copying mechanism
width_view = np.pi/3 # Width of the field of view
SEED = 2302 # Seed for the random number generation


# We create a bird list
birdlist = [viv.BirdVision(width_view) for _ in range(N_birds)]

# We set the coordinates of the birds to random
for i in range(N_birds):
    birdlist[i].x = length*np.random.rand() # Random x coordinate in [0, length)
    birdlist[i].y = length*np.random.rand() # Random y coordinate in [0, length)
    birdlist[i].theta = 2*np.pi*(np.random.rand() - 0.5) # Random angle in [-pi, pi)

# We create the simulation object
sim = viv.VicsekSimulation(birdlist, L=length, v=vel, eta=ang_noise, seed=SEED)

# We save the arrays with the coordinates. The arrays don't own the data
# they just acces the data that's stored in C++. As a results these arrays
# will be updated as the simulation goes on
x, y, theta = sim.x, sim.y, sim.theta


# We configure the figure to create the GIF
fig, ax = plt.subplots()
fig.set_tight_layout(True)

ax.set_xlim(0, length)
ax.set_ylim(0, length)

ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)

quivBirds = ax.quiver(x, y, np.cos(theta), np.sin(theta), color="#101062")

ax.set_title("t = 0")

# We define the function that updates to the next frame of the GIF
# the parameter n is the number of the current frame
def update(n: int) -> None:
    ax.set_title(f"t = {n*N_steps}")
    sim.update(N_steps) # Updates the simulaton
    quivBirds.set_UVC(np.cos(theta), np.sin(theta))
    quivBirds.set_offsets(np.c_[x, y])

# We create the GIF
anim = FuncAnimation(fig, update, frames=N_frames, interval=100)
anim.save("test_anim.gif", dpi=80, writer="imagemagick")
