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
pyvicsek.Bird(self: pyvicsek.Bird, x: float = 0, y: float = 0, theta: float = 0)
```

Sets a bird in the position $\left(x, y\right)$ and with angle theta.

**Way 2 of constructing a Bird**

```Python
pyvicsek.Bird(self: pyvicsek.Bird, sim: pyvicsek.VicsekSimulation)
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

### pyvicsek.VicsekSimulation

This class is the key one for the simulations. It contains all the information needed to run a simulation. Internally it stores a list of Bird objects that are used to run the simulation. It also stores all the simulation parameters. In order to optimize the simulation a link list is created, this list plus some information about the neighbours of each cell in the link list is also stored in this class.

#### Attribitues

All the attributes except for the simulation time are readonly since changing them will change the simulation completely. However if you really want to change the parameters you can use the setters. The accesible attributes of this class are:

1. **t**: Time of the simulation. It is automatically increased when we call the method ``VicsekSimulation.update``. It does not affect the simulation and can be set freely.
2. **N**: Number of birds of the simulation. It can be modified using ``VicsekSimulation.set_N`` and ``VicsekSimulation.push``.
3. **L**: Length of the box in which the simulation takes place. For some simulations the density of birds may be an important factor. Even though the density is not an attribute and cannot be set keep in mind that it is just $\rho = N/L^2$. The length can only be modified using ``VicsekSimulation.set_L``.
4. **exact_cells**: Boolean variable indicating if the length is an exact number of the radius of vision of the bird (so basically if the length is an integer). It is modified when using ``VicsekSimulation.set_L``. The simulation will run faster if this attribute is ``True``, so it is always nice to check.
5. **eta**: Angular noise of the copying mechanism. It can only be modified via ``VicsekSimulation.set_eta``.
6. **v**: Speed of the birds, it can be modified using ``VicsekSimulation.set_v``.
7. **x**, **y**, **theta**: These are numpy arrays containing all the coordinates and angles of the Birds in the simulation. They are readonly and the data is stored in C++ and accesed directly from Python. As a result no copy of the data from C++ to Python is required.

#### Constructors

There are two possible constructors for this class.

**Constructor that receives the number of birds**:

```Python
pyvicsek.VicsekSimulation(self: pyvicsek.VicsekSimulation, N: int, v: float = 0.03, eta: float = 2.0, L: float = 7.0, seed: int = 0)
```

N is the number of birds that we want to simulate. These birds will be created from the Bird base class with random coordinates. Then v is the speed of the birds, eta is the angluar noise in the copying mechanism, L is the box length and seed is a seed use for the random number generation. By defualt seed is set to 0 so if you run the exact same simulation twice the results will be the same.

**Constructor that receives a list of birds:**

```Python
pyvicsek.VicsekSimulation(self: pyvicsek.VicsekSimulation, birds: List[pyvicsek.Bird], v: float = 0.03, eta: float = 2.0, L: float = 7.0, seed: int = 0)
```

All the optional parameters are the same ones but now instead of giving a number of birds we have to give a list with all the birds that can also be objects derived from the base class Bird. The parameters of the birds will be the ones from the last simulation in which they are included. This means that it is allow to first add a bird to a simulation sim1 and then add it to the simulation sim2. But then trying to update sim1 can lead to wrong results and even crash the simulation.

#### Methods

The pyvicsek.VicsekSimulation class has the following methods:

```Python
pyvicsek.VicsekSimulation.update(self: pyvicsek.VicsekSimulation, N_steps: int = 1)
```

Updates the position of the birds N_steps times. The time attribute will be increased by N_steps. The data about the simulation is not saved until the last step is finished, as a result, it is more efficient to call this method one time with the desired number of steps instead of calling it several times. The general strategy should be: simulation.update(SomeSteps) :arrow_right: Do measures :arrow_right: simulation.update(SomeSteps) :arrow_right: Do measures again :arrow_right: ...

```Python
pyvicsek.VicsekSimulation.set_v(self: pyvicsek.VicsekSimulation, v: float)
```

This method modifies the v attribute and, as a result, sets the speed of all the birds in the simulation

```Python
pyvicsek.VicsekSimulation.set_eta(self: pyvicsek.VicsekSimulation, eta: float)
```

This method modifies the eta attribute thus changing the angular noise of the copying mechanism.

```Python
pyvicsek.VicsekSimulation.set_L(self: pyvicsek.VicsekSimulation, L: float)
```

This method sets the length of the box in which the simulation takes place. If some birds are as a result out of bounds they will be moved to fit into the new box. This method also modifies the exact_cells attribute. As mentioned before this attribute will be true if L is an integer and in that case the simulation will run faster.

```Python
pyvicsek.VicsekSimulation.set_N(self: pyvicsek.VicsekSimulation, N: int)
```

Sets the number of birds of the simulation. If the new number is bigger more objects of the Bird class will be created with random coordinates, otherwise the list of birds in the simulation will just shrink.

```Python
pyvicsek.VicsekSimulation.push(self: pyvicek.VicsekSimulation, birds: pyvicsek.Bird)
pyvicsek.VicsekSimulation.push(self: pyvicek.VicsekSimulation, birds: List[pyvicsek.Bird])
```

This is an overloaded method. It allows you to add a single bird (or derived object) to the simulation or a whole list of birds. The number of birds will grow as a consequence. As it has been already mentioned the same bird object should not belong to more than one simulation in progress. If this happens the results will be wrong and the simulation can even crash.

```Python
pyvicsek.VicsekSimulation.randomize(self: pyvicek.VicsekSimulation)
```

Set all the birds to random positions and set all their directions to random angles.

```Python
pyvicsek.VicsekSimulation.compute_va(self: pyvicek.VicsekSimulation) -> float
```

Computes the order parameter of the simulation. The definition of this parameter has already been explained in the part about the Vicsek model, but this parameter will be 0 if all the directions are completely random and 1 if all the directions are the same.

### Examples

Some examples are created in the ``examples`` directory. The first and most basic one of them is just for obtaining the order parameter:

```Python
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
```

## Defining new bird behaviour

This module is created in such a way that adding a new type of bird with different behaviour is more or less easy. A full example is shown in the directory *derived_bird_example*. Let's explain a little bit this example.

The goal is to create a new kind of bird that can only copy other birds if they are at a distance smaller than its vision radius ($r=1$) **and** they are inside of the field of view. The width of the field of view will be a new parameter that we need to specify.

In order to create this new bird and bring it to live in Python we need to basically perform two steps. First of all we will create a new C++ class that inherits from the Bird class. The header for this new class could be defined as follows:

```C++
#pragma once

#include <cmath>
#include <list>
#include <numbers>
#include <memory>

#include "vicsek.hpp"

// We create a new class the inherits from Bird
class BirdVision : public Bird {
	double ang_field_view; // Field of view angle in radians
protected:
	void sum_birds_in_view(const std::list< std::shared_ptr<Bird> >&, double&, double&) const;
public:
	BirdVision(double FieldOfView = std::numbers::pi/3, double x = 0, double y = 0, double theta = 0)
	: Bird(x, y, theta), ang_field_view(FieldOfView) {}
	BirdVision(VicsekSimulation* sim, double FieldOfView = std::numbers::pi/3)
	: Bird(sim), ang_field_view(FieldOfView) {}

	virtual ~BirdVision() = default;

	// Gets the angle between this bird and the bird given by the pointer
	double get_angle_to(const std::shared_ptr<Bird>&) const;

	// Overwrite the average angle function
	virtual double average_angle(const std::vector< std::vector<unsigned> >&, const std::vector< std::list< std::shared_ptr<Bird> > >&) const;
};
```
