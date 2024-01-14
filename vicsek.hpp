#pragma once

#include <vector>
#include <list>
#include <Eigen/Dense>
#include <memory>
#include <algorithm>
#include <random>
#include <cmath>
#include <numbers>
#include <ostream>

class VicsekSimulation;

class Bird {
	// Position
	double x, y;
	// Angle
	double theta;

	// Link cell indeces
	unsigned i, j;

	// Variables to do a partial update
	double x_new, y_new, theta_new;

	// Pointer to the simulation it belongs
	VicsekSimulation* sim;

protected:
	void sum_birds_in(const std::list< std::shared_ptr<Bird> >&, double&, double&) const;

public:
	// Default constructor
	Bird(VicsekSimulation*);
	Bird(double, double, double);

	// Virtual default destructor
	virtual ~Bird() = default;

	// Copy constructor
	Bird(const Bird& o)
	: x(o.x), y(o.y), theta(o.theta), i(o.i), j(o.j), x_new(o.x_new), y_new(o.y_new), theta_new(o.theta_new), sim(o.sim) {}

	// Setters
	void set_x(double);
	void set_y(double);
	void set_pos(double, double);
	void set_pos_new(double, double);
	void set_angle(double angle) {
		theta = std::fmod(angle, 2*std::numbers::pi);

		if (theta < 0) theta += 2*std::numbers::pi;
	}
	void set_angle_new(double angle) {
		theta_new = std::fmod(angle, 2*std::numbers::pi);

		if (theta_new < 0) theta_new += 2*std::numbers::pi;
	}
	void set_sim(VicsekSimulation* s) {
		sim = s;
		set_pos(x, y);
	}

	void randomize();


	// Getters
	double get_x() const { return x; }
	double get_y() const { return y; }
	double get_angle() const { return theta; }
	unsigned get_i() const { return i; }
	unsigned get_j() const { return j; }
	VicsekSimulation* get_sim_ptr() const { return sim; }

	// Compute the square of the distance using periodic boundary conditions
	double distance_sq(const Bird&) const;

	// Computes the average angle fo the birds in the list
	// that are at a distance r from this bird
	virtual double average_angle(const std::vector< std::vector<unsigned> >&, const std::vector< std::list< std::shared_ptr<Bird> > >&) const;

	// Updates but instead of updating it directly
	// it saves everything in a new variable
	virtual void update_to_new(const std::vector< std::vector<unsigned> >&, const std::vector< std::list< std::shared_ptr<Bird> > >&);

	// It loads the content of the new variables
	// into the usual variables
	void update_from_new(std::vector< std::list< std::shared_ptr<Bird> > >&);
};

class VicsekSimulation {
	// Attributes
	unsigned N_birds; // Number of birds
	double vel; // Modulus of the velocity of the birds
	double eta; // Noise introduced in the copying angle mechanism
	std::list< std::shared_ptr<Bird> > bird_list; // List containing the birds
	double L; // Length of the box
	unsigned N_cells; // Number of cells in each direction
	int time;

	bool exact_cells;
	std::vector< std::vector<unsigned> > neighbours;
	std::vector< std::list< std::shared_ptr<Bird> > > bird_linklist;

	// Attributes for random number generation
	std::mt19937 gen;
	std::uniform_real_distribution<double> ran_theta;
	std::uniform_real_distribution<double> ran_u;

	Eigen::VectorXd x, y, theta;

	// Private methods
	void compute_neighbours();
	void compute_linklist();
public:
	// Default constructor
	VicsekSimulation(unsigned, double, double, double, int);
	VicsekSimulation(const std::list< std::shared_ptr<Bird> >&, double, double, double, int);

	// Getters
	unsigned get_N() const { return N_birds; }
	unsigned get_N_cells() const { return N_cells; }
	double get_v() const { return vel; }
	double get_eta() const { return eta; }
	double get_L() const { return L; }
	bool get_exact_cells() const { return exact_cells; }
	double get_ran_u() { return ran_u(gen); }
	double get_ran_theta() { return ran_theta(gen); }
	int get_time() const { return time; }
	const Eigen::VectorXd& get_x() const { return x; }
	const Eigen::VectorXd& get_y() const { return y; }
	const Eigen::VectorXd& get_angle() const { return theta; }

	// Setters
	void set_time(int t) { time = t; }
	void set_N(unsigned);
	void set_v(double v) { vel = v; }
	void set_eta(double noise) {
		eta = noise;
		ran_theta = std::uniform_real_distribution<double>(-eta/2, eta/2);
	}
	void set_L(double);

	void push(const std::shared_ptr<Bird>&);
	void push(const std::list< std::shared_ptr<Bird> >&);

	void randomize();

	void change_bird_in_linklist(Bird*, unsigned, unsigned);

	// Get the position in a 1D vector given the indicies (i, j)
	unsigned pos_by_indeces(unsigned i, unsigned j) const { return i*N_cells + j; }

	// Update the position of the birds a given number of steps (by default one)
	void update(unsigned);

	// Computes the order parameter (v_a): 1 -> all angles are the same and 0 -> completely random angles. This parameter is only computed for the current state, typically one would like to compute the time average
	double compute_va() const;

	// Operator overloading to show the current state of the simulation
	friend std::ostream& operator<<(std::ostream&, const VicsekSimulation&);
};
