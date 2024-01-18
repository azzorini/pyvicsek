#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>
#include <list>
#include <memory>

#include "vicsek.hpp"
#include "vicsek_vision.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pyvicsekvision, m) {
	m.doc() = "pyvicsek is a Python module to run Vicsek simulations efficiently";

	py::class_< Bird, std::shared_ptr<Bird> >(m, "Bird")
		.def(py::init<double, double, double>(), "Creates a bird with the given coordinates (x, y, θ). Where (x, y) are the cartesian coordinates and θ is the bird's angle.", py::arg("x") = 0, py::arg("y") = 0, py::arg("theta") = 0)
		.def(py::init<VicsekSimulation*>(), "Creates a bird with a random angle and coordinates randomly distributed in the interval [0, L), where L is the length of the VicsekSimulation passed to the constructor.\n\nWarning: The constructor does not add the Bird to the simulation if you want to do that you will need to use the method VicsekSimulation.push.", py::arg("sim"))
		.def_property("x", &Bird::get_x, &Bird::set_x, "Cartesian coordinate of the bird in the horizontal direction.")
		.def_property("y", &Bird::get_y, &Bird::set_y, "Cartesian coordinate of the bird in the vertical direction.")
		.def_property("theta", &Bird::get_angle, &Bird::set_angle, "Angle of the bird.")
		.def_property_readonly("sim", &Bird::get_sim_ptr, "Object that contains all the simulation parameters of the bird. The bird coordinates will be forcefully at the interval [0, L].")
		.def("__repr__",
		     [](const Bird& b) {
			     return '(' + std::to_string(b.get_x()) + ", " + std::to_string(b.get_y()) + ", " + std::to_string(b.get_angle()) + ')';
		     }, "String representation of a bird: (x, y, θ).");

	// We add the bindings for the new Bird type
	py::class_< BirdVision, Bird, std::shared_ptr<BirdVision> >(m, "BirdVision")
		.def(py::init<double, double, double, double>(), "Creates a bird with a field of view width given by WidthView and the given coordinates (x, y, theta). Where (x, y) are the cartesian coordinates and theta is the bird's angle.", py::arg("WidthView") = std::numbers::pi/3, py::arg("x") = 0, py::arg("y") = 0, py::arg("theta") = 0)
		.def(py::init<VicsekSimulation*, double>(), "Creates a bird with a field of view width given by WidthView, random angle and random coordinates randomly distributed in the interval [0, L), where L is the length of the VicsekSimulation passed to the constructor.\n\nWarning: The constructor does not add the Bird to the simulation if you want to do that you will need to use the method VicsekSimulation.push.", py::arg("sim"), py::arg("WidthView") = std::numbers::pi/3)
		.def_property("x", &BirdVision::get_x, &BirdVision::set_x, "Cartesian coordinate of the bird in the horizontal direction.")
		.def_property("y", &BirdVision::get_y, &BirdVision::set_y, "Cartesian coordinate of the bird in the vertical direction.")
		.def_property("theta", &BirdVision::get_angle, &BirdVision::set_angle, "Angle of the bird.")
		.def_property_readonly("sim", &BirdVision::get_sim_ptr, "Object that contains all the simulation parameters of the bird. The bird coordinates will be forcefully at the interval [0, L].")
		.def("__repr__",
		     [](const BirdVision& b) {
			     return '(' + std::to_string(b.get_x()) + ", " + std::to_string(b.get_y()) + ", " + std::to_string(b.get_angle()) + ')';
		     }, "String representation of a bird: (x, y, θ).");

	py::class_<VicsekSimulation>(m, "VicsekSimulation")
		.def(py::init<unsigned, double, double, double, int>(), "Builds a simulation with the following parameters:\n\n\t-N: Numbers of Birds in the simulation.\n\t-v: Speed of the birds. The time units increased by one in every update and the length units are chosen so the radius of vision of the bird is 1. Default value: v = 0.03.\n\t-eta: Angular noise in the copying mechanism. Every time the angle is updated a noise between (-eta/2, eta/2) will be introduced. Default value: eta = 2.0.\n\t-L: Length of the box in which the simulation takes place. The radius of vision of the bird is 1 length unit. Periodic boundary conditions are assumed. Warning: if the length is a integer the simulation will be faster, even though fractional lengths are also allowed. Default value: L = 7.\n\t-seed: Seed for the random number generator. Default value: seed = 0.",  py::arg("N"), py::arg("v") = .03, py::arg("eta") = 2., py::arg("L") = 7., py::arg("seed") = 0)
		.def(py::init<const std::list< std::shared_ptr<Bird> >&, double, double, double, int>(), "All the optional arguments are the same as the ones explained for the other constructor but now the only required argument birds is a list of Bird objects (objects derived from Bird are also allowed)", py::arg("birds"), py::arg("v") = .03, py::arg("eta") = 2., py::arg("L") = 7., py::arg("seed") = 0)
		.def_property("t", &VicsekSimulation::get_time, &VicsekSimulation::set_time, "Simulation time. It does not have an effect on the simulation and can be modified.")
		.def_property_readonly("x", &VicsekSimulation::get_x, "Array with all the x coordinates of the birds", py::return_value_policy::reference_internal)
		.def_property_readonly("y", &VicsekSimulation::get_y, "Array with all the y coordinates of the birds", py::return_value_policy::reference_internal)
		.def_property_readonly("theta", &VicsekSimulation::get_angle, "Array with all the angles of the birds", py::return_value_policy::reference_internal)
		.def_property_readonly("L", &VicsekSimulation::get_L, "Length of the box in which the simulation takes place. It can only be modified using VicsekSimulation.set_L.")
		.def_property_readonly("v", &VicsekSimulation::get_v, "Speed of the birds. It can only be modified using VicsekSimulation.set_v.")
		.def_property_readonly("N", &VicsekSimulation::get_N, "Number of birds in the simulation. It can only be modified using VicsekSimulation.set_N.")
		.def_property_readonly("eta", &VicsekSimulation::get_eta, "Angular noise in the copying mechanism. It can only be modified using VicsekSimulation.set_eta.")
		.def_property_readonly("exact_cells", &VicsekSimulation::get_exact_cells, "Boolean variable indicating if the length is an exact number of the radius of vision of the bird (so basically if the length is an integer). It is modified when using VicsekSimulation.set_L. The simulation will run faster if this attribute is True.")
		.def("set_N", &VicsekSimulation::set_N, "Sets the number of birds in the simulation. If the new number of birds is bigger new Bird objects will be generated.", py::arg("N"))
		.def("set_L", &VicsekSimulation::set_L, "Sets the length of the box in which the simulation takes place.", py::arg("L"))
		.def("set_eta", &VicsekSimulation::set_eta, "Sets the angular noise of the copying mechanism.", py::arg("eta"))
		.def("set_v", &VicsekSimulation::set_v, "Sets the speed of the birds.", py::arg("v"))
		.def("push", py::overload_cast< const std::shared_ptr<Bird>& >(&VicsekSimulation::push), "Adds a single Bird (or derived object) to the simulation.",  py::arg("birds"))
		.def("push", py::overload_cast< const std::list< std::shared_ptr<Bird> >& >(&VicsekSimulation::push), "Adds a list of birds (or derived objects) to the simulation.", py::arg("birds"))
		.def("randomize", &VicsekSimulation::randomize, "Sets all the coordinates of the birds to random values.")
		.def("update", &VicsekSimulation::update, "Updates the simulation N_steps times and increases the time consequently. If no measures are going to be taken it is more efficient to give a higher value for N_steps.", py::arg("N_steps") = 1)
		.def("compute_va", &VicsekSimulation::compute_va, "Computes the order parameter which is just the modulus of the sum over birds of exp(iθ_j) where θ_j is the angle of each bird.")
		.def("__repr__",
		     [](const VicsekSimulation& sim) {
		     return "<pyvicsek.VicsekSimulation with N=" + std::to_string(sim.get_N()) + ", L=" + std::to_string(sim.get_L()) + ", eta=" + std::to_string(sim.get_eta()) + '>';
		     }, "Representation of the simulation.");
}
