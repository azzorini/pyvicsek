# pyvicsek

A Python module for efficient Vicsek simulations.

## About it

This is a Python module made in C++ using [pybind11](https://github.com/pybind/pybind11). The standard used to compile the module is C++20. A Makefile is given to compile the module the only two requisites are:

1. Python needs to be installed in the computer. ``python`` and ``python-config`` have to be in the PATH.
2. Pybind11 needs to be installed.

## Vicsek Simulation

This purpose of this module is to be able to run a simulation of the Vicsek model. So let's start by explainig what is the Vicsek model. This model tries to replicate flocks of birds.

So in the model we will have $N$ birds. These birds will move in a box of length $L$. For the simulation periodic boundary conditions are used. Each bird has a radius of vision $r$. Without any loss of generality we set $r=1$. This means that the unit of length is just the radius of vision. All the birds are given the same constant speed $v$. At every time step every single bird will move a distance $v\Delta t$ (where $\Delta t = 1$ fixing now the time units) in whatever direction it is pointing. Then the direction of the bird is also updated. In order to update its direction the bird will try to imitate the average direction of all the birds included in its radius of vision (including itself). But to this average we will add a noise which is just a random number between $\left(-\eta / 2,\ \eta / 2\right)$. Where $\eta$ is another parameter of the model.

A sketch of how a simulation looks is shown in the following image:

![Vicsek model sketch](https://github.com/azzorini/pyvicsek/blob/main/img/VicsekModelSketch.png?raw=true)

We will be able to define and modify all the explained parameters.

Now let's briefly described the mathematical part of the model. At every time step we update the bird $i$ using the following equation:

$$
\boldsymbol{x}_i(t+1) = \boldsymbol{x}_i(t) + \boldsymbol{v}_i(t)\Delta t
$$

The velocity $\boldsymbol{v}_i$ of the bird is:

$$
\boldsymbol{v}_i = v\left[ \hat{\boldsymbol{x}}\cos\theta_i(t) + \hat{\boldsymbol{y}}\sin\theta_i(t)\right]
$$

Then the angle is also updated as follows:

$$
\theta_i(t+1) = \left\langle \theta(t) \right\rangle_r + \Delta\theta
$$

where $\Delta\theta$ is a random number in the interval $\left(-\eta / 2,\ \eta / 2\right)$ and $\left\langle \cdots \right\rangle_r$ denotes the average over the birds that are at a distance smaller than $r$. The best way to do the average is not so trivial. We are doing the averga in the following way:

$$
\left\langle \theta(t) \right\rangle_r = \arctan \left[ \frac{\left\langle \sin\theta(t)\right\rangle_r}{\left\langle \cos\theta(t) \right\rangle_r}\right]
$$

The advantage of taking the average in this way is that it will work as expected even if we have angles that are not in the interval $\left[-\pi,\ \pi\right]$. This interval is chosen since the version of the arctangent that we are using (``std::atan2`` from ``<cmath>``) gives us results in that particular interval.

To see if all the birds are pointing to the same direction or if they are pointing random directions we will use an order parameter. This parameter will be $0$ if all the birds are randomly directed and it will be $1$ if they are all in the same direction. The order parameter is defined as

$$
v_a(t) = \frac{1}{N}\left| \sum_{j=1}^N e^{i\theta_j(t)} \right|
$$

The model now has already been explained let's jump to the coding part.

## Python API

There are two classes that we need to understand in order to work with this module: pyvicsek.Bird and pyvicsek.VicsekSimulation.

### pyvicsek.Bird

This class represents a single bird. It contains information about the bird's coordinates and also about the simulation that sets the bird parameters.

The attributes of the Bird class are:

1. **x**: Cartesian coordinate of the bird in the horizontal direction.
2. **y**: Cartesian coordinate of the bird in the vertical direction.
3. $\boldsymbol{\theta}$: Angle of the bird.
4. **sim**: pyvicsek.VicsekSimulation object that contains all the simulation parameters. This attribute cannot be directly modified but when the bird is added to a simulation this parameter is changed. The coordinates of the birds are restricted to the interval $\left[0,\ L\right]$ where $L$ is defined in the VicsekSimulation object.

There are two ways of creating a bird object:

**Way 1 of constructing a Bird**

```Python
pyvicsek.Bird(x: float = 0, y: float = 0, theta: float = 0)
```

Sets a bird in the position $\left(x, y\right)$ and with angle theta.

**Way 2 of constructing a Bird**

```Python
pyvicsek.Bird(sim: pyvicsek.VicsekSimulation)
```

Sets the attribute sim of the Bird, sets the coordinate at some random value in the interval $\left[0, L\right]$ where $L$ is taken from the sim parameter and sets the angle at a random value.

The representation of the bird is given by: (x, y, theta). Example:

```Python
import pyvicsek as vi

b = vi.Bird(2, 3, 1.5)

print(b) # -> It will show (2.000000, 3.000000, 1.500000)

b.x += 2.5

print(b) # -> It will show (4.500000, 3.000000, 1.500000)
```
