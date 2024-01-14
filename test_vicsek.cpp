#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <list>
#include <memory>
#include <numbers>

#include "vicsek.hpp"

const int SEED = 2802;

int main(int argc, char* argv[]) {
	unsigned i;

	if (argc < 4) {
		std::cout << "You need to put at least three arguments: Number of birds, length of the box and noise strength\n";
		return 1;
	}

	const unsigned N_birds = std::stoi(argv[1]);
	const double L = std::stod(argv[2]);
	const double eta = std::stod(argv[3]);

	const unsigned N = 400;
	const unsigned N_steps = 5;
	const double vel = .03;

	std::mt19937 gen(SEED);
	std::uniform_real_distribution<double> ran_u(0, 1);

	std::list< std::shared_ptr<Bird> > birds;

	std::cout << "Running simulation with N_birds = " << N_birds << ", L = " << L << " and eta = " << eta << '\n';

	std::ofstream fbirds("birds_" + std::to_string(N_birds) + "_N_" + std::to_string(N) + "_L_" + std::to_string(L) + "_test.txt");
	std::ofstream fva("va_" + std::to_string(N_birds) + "_N_" + std::to_string(N) + "_L_" + std::to_string(L) + "_test.txt");

	fbirds << "#t\tx\ty\ttheta\n";
	fva << "#t\tva\n";

	// VicsekSimulation sim(N_birds, vel, eta, L, SEED); // First constructor

	// Second constructor
	for (i = 0; i < N_birds; i++) {
		birds.push_back(std::make_shared<Bird>(L*ran_u(gen), L*ran_u(gen), 2*std::numbers::pi*ran_u(gen)));
	}
	VicsekSimulation sim(birds, vel, eta, L, SEED+3*N_birds);

	// Simulation loop
	for (i = 0; i < N; i++) {
		sim.update(N_steps);

		fbirds << sim << '\n';

		fva << sim.get_time() << '\t' << sim.compute_va() << '\n';
	}

	fbirds.close();
	fva.close();
}
