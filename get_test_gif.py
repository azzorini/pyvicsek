#!/usr/bin/python

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def update_arrays(f, x, y, angle):
    line = f.readline()
    i = 0
    while line != "\n" and line != "":
        if line[0] != '#':
            x[i], y[i], angle[i] = (float(s) for s in line.split()[1:])
            i += 1
        line = f.readline()

def update(n):
    ax.set_title(f"t = {n*N_steps}")
    update_arrays(f, x, y, theta)
    q.set_UVC(np.cos(theta), np.sin(theta))
    q.set_offsets(np.c_[x, y])

for filename in os.listdir():
    if ".txt" in filename and "birds_" in filename and "_test" in filename:
        print(f"Creating gif for {filename}...")
        name_split = filename.split('_')
        N, t_f, L = int(name_split[1]), int(name_split[3]), float(name_split[5])
        N_steps = 5
        fig, ax = plt.subplots()
        fig.set_tight_layout(True)

        ax.set_xlim(0, L)
        ax.set_ylim(0, L)

        ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False,
                        right=False, left=False, labelleft=False)

        x, y, theta = .5*L*np.ones(N), .5*L*np.ones(N), np.zeros(N)
        f = open(filename, "r")
        q = ax.quiver(x, y, np.cos(theta), np.sin(theta), color="#101062")

        ax.set_title("t = 0")

        anim = FuncAnimation(fig, update, frames=t_f, interval=100)
        anim.save(filename.replace(".txt", ".gif"), dpi=80, writer="imagemagick")
