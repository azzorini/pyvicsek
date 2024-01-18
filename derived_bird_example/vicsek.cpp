#include "vicsek.hpp"

// CLASS: Bird

void Bird::sum_birds_in(const std::list< std::shared_ptr<Bird> >& birds, double& sum_sin, double& sum_cos) const {
	for (const std::shared_ptr<Bird>&  p : birds) {
		if (distance_sq(p) < 1) {
			sum_sin += std::sin(p->theta);
			sum_cos += std::cos(p->theta);
		}
	}
}

Bird::Bird(double x_, double y_, double ang)
: x(x_), y(y_), theta(ang), x_new(0), y_new(0), theta_new(0), sim(nullptr) {
	i = (unsigned) x;
	j = (unsigned) y;
}

Bird::Bird(VicsekSimulation* sim_ptr)
: x_new(0), y_new(0), theta_new(0), sim(sim_ptr) {
	const double L = sim->get_L();

	theta = 2*std::numbers::pi*(sim->get_ran_u() - 0.5);

	x = L*sim->get_ran_u();
	y = L*sim->get_ran_u();

	i = (unsigned) x;
	j = (unsigned) y;
}

void Bird::set_x(double x_) {
	double L;
	unsigned old_i;

	x = x_;
	old_i = i;

	if (sim != nullptr) {
		L = sim->get_L();

		while (x < 0) x += L;
		while (x > L) x -= L;
		i = (unsigned) x;

		if (old_i != i) {
			sim->change_bird_in_linklist(this, old_i, j);
		}

		return;
	}

	i = (unsigned) x;
}

void Bird::set_y(double y_) {
	double L;
	unsigned old_j;

	y = y_;

	old_j = j;

	if (sim != nullptr) {
		L = sim->get_L();

		while (y < 0) y += L;
		while (y > L) y -= L;
		j = (unsigned) y;

		if (old_j != j) {
			sim->change_bird_in_linklist(this, i, old_j);
		}
	}

	j = (unsigned) y;
}

void Bird::set_pos(double x_, double y_) {
	double L;

	x = x_;
	y = y_;

	if (sim != nullptr) {
		L = sim->get_L();
		while (x < 0) x += L;
		while (x > L) x -= L;
		while (y < 0) y += L;
		while (y > L) y -= L;
	}

	i = (unsigned) x;
	j = (unsigned) y;
}

void Bird::set_pos_new(double x_, double y_) {
	const double L = sim->get_L();

	x_new = x_;
	y_new = y_;

	while (x_new < 0) x_new += L;
	while (x_new > L) x_new -= L;
	while (y_new < 0) y_new += L;
	while (y_new > L) y_new -= L;
}

void Bird::randomize() {
	const double L = sim->get_L();

	x = L*sim->get_ran_u();
	y = L*sim->get_ran_u();
	theta = 2*std::numbers::pi*(sim->get_ran_u() - 0.5);

	i = (unsigned) x;
	j = (unsigned) y;
}

double Bird::distance_sq(const std::shared_ptr<Bird>& p) const {
	const double L = sim->get_L();

	double dx = std::abs(x - p->x);
	double dy = std::abs(y - p->y);

	if (dx > .5*L) dx = L - dx;
	if (dy > .5*L) dy = L - dy;

	return dx*dx + dy*dy;
}

double Bird::average_angle(const std::vector< std::vector<unsigned> >& neighbours, const std::vector< std::list< std::shared_ptr<Bird> > >& link_list) const {
	unsigned ind = sim->pos_by_indeces(i, j);
	double sum_sin = .0, sum_cos = .0;
	unsigned n, N_neigh = neighbours[0].size(), N_cells;

	sum_birds_in(link_list[ind], sum_sin, sum_cos);

	for (n = 0; n < N_neigh; n++) {
		sum_birds_in(link_list[neighbours[ind][n]], sum_sin, sum_cos);
	}

	// Some additional checks if the number of cells does not exactly match the length
	if (!sim->get_exact_cells()) {
		N_cells = sim->get_N_cells();
		if (i == 0) {
			n = sim->pos_by_indeces(N_cells-1, j);
			sum_birds_in(link_list[neighbours[n][1]], sum_sin, sum_cos);
			sum_birds_in(link_list[neighbours[n][4]], sum_sin, sum_cos);
			sum_birds_in(link_list[neighbours[n][5]], sum_sin, sum_cos);
		} else if (i == N_cells - 2) {
				n = sim->pos_by_indeces(N_cells-1, j);
				sum_birds_in(link_list[neighbours[n][3]], sum_sin, sum_cos);
				sum_birds_in(link_list[neighbours[n][6]], sum_sin, sum_cos);
				sum_birds_in(link_list[neighbours[n][7]], sum_sin, sum_cos);

		}
		if (j == 0) {
			n = sim->pos_by_indeces(i, N_cells-1);
			sum_birds_in(link_list[neighbours[n][0]], sum_sin, sum_cos);
			sum_birds_in(link_list[neighbours[n][4]], sum_sin, sum_cos);
			sum_birds_in(link_list[neighbours[n][7]], sum_sin, sum_cos);
		} else if (j == N_cells - 2) {
				n = sim->pos_by_indeces(i, N_cells-1);
				sum_birds_in(link_list[neighbours[n][2]], sum_sin, sum_cos);
				sum_birds_in(link_list[neighbours[n][5]], sum_sin, sum_cos);
				sum_birds_in(link_list[neighbours[n][6]], sum_sin, sum_cos);

		}

	}

	// Technically what we should do is the average of the sums
	// but since std::atan2 only care about the proportion
	// there is no point in doing the averages
	return std::atan2(sum_sin, sum_cos);
}

void Bird::update_to_new(const std::vector< std::vector<unsigned> >& neighbours, const std::vector< std::list< std::shared_ptr<Bird> > >& link_list) {
	const double L = sim->get_L();
	const double v = sim->get_v();

	// We update the position
	x_new = x + v*std::cos(theta);
	y_new = y + v*std::sin(theta);

	// We check that the bird has not crossed the borders
	// otherwise we put it in the right place
	if (x_new < 0) x_new += L;
	if (x_new > L) x_new -= L;

	if (y_new < 0) y_new += L;
	if (y_new > L) y_new -= L;

	// We update the angle
	theta_new = average_angle(neighbours, link_list) + sim->get_ran_theta();
}

void Bird::update_from_new(std::vector< std::list< std::shared_ptr<Bird> > >& link_list) {
	unsigned aux, ind = sim->pos_by_indeces(i, j);
	bool changed = false;
	std::list< std::shared_ptr<Bird> >::iterator it;

	x = x_new;
	y = y_new;
	theta = theta_new;

	// We check if we need to update i and j
	aux = (unsigned) x;
	if (i != aux) {
		i = aux;
		changed = true;
	}

	aux = (unsigned) y;
	if (j != aux) {
		j = aux;
		changed = true;
	}

	if (changed) {
		// We update the link_list

		// We assume that this list is ok and we did find it
		it = std::find_if(link_list[ind].begin(), link_list[ind].end(), [&](const std::shared_ptr<Bird>& ptr){ return ptr.get() == this; });
		// We add the pointer to the new position
		link_list[sim->pos_by_indeces(i, j)].push_back(*it);
		// We delete de pointer from the old position
		link_list[ind].erase(it);
	}
}

// CLASS: VicsekSimulation

void VicsekSimulation::compute_neighbours() {
	unsigned i, j, pos;
	unsigned ant_i, next_i, ant_j, next_j;

	for (i = 0; i < N_cells; i++) {
		ant_i = (i > 0 ? i-1 : N_cells-1);
		next_i = (i < N_cells-1 ? i+1 : 0);
		for (j = 0; j < N_cells; j++) {
			ant_j = (j > 0 ? j-1 : N_cells-1);
			next_j = (j < N_cells-1 ? j+1 : 0);
			pos = pos_by_indeces(i, j);
			neighbours[pos][0] = pos_by_indeces(ant_i, j);
			neighbours[pos][1] = pos_by_indeces(i, ant_j);
			neighbours[pos][2] = pos_by_indeces(next_i, j);
			neighbours[pos][3] = pos_by_indeces(i, next_j);
			neighbours[pos][4] = pos_by_indeces(ant_i, ant_j);
			neighbours[pos][5] = pos_by_indeces(next_i, ant_j);
			neighbours[pos][6] = pos_by_indeces(next_i, next_j);
			neighbours[pos][7] = pos_by_indeces(ant_i, next_j);
		}
	}
}

void VicsekSimulation::compute_linklist() {
	for (const std::shared_ptr<Bird>& ptr : bird_list)
		bird_linklist[pos_by_indeces(ptr->get_i(), ptr->get_j())].push_back(ptr);
}

VicsekSimulation::VicsekSimulation(unsigned BirdNumber, double velocity, double ang_noise, double length, int seed)
: N_birds(BirdNumber), vel(velocity), eta(ang_noise), L(length), N_cells(std::ceil(length)), time(0) {
	unsigned i, N_cells_sq = N_cells*N_cells;
	std::shared_ptr<Bird> ptr;

	exact_cells = !(std::ceil(length) > length);
	gen = std::mt19937(seed);
	ran_theta = std::uniform_real_distribution<double>(-eta/2, eta/2);
	ran_u = std::uniform_real_distribution<double>(0, 1);
	x = Eigen::VectorXd(N_birds);
	y = Eigen::VectorXd(N_birds);
	theta = Eigen::VectorXd(N_birds);

	for (i = 0; i < N_cells_sq; i++) {
		neighbours.push_back(std::vector<unsigned>(8));
		bird_linklist.push_back(std::list< std::shared_ptr<Bird> >());
	}

	for (i = 0; i < N_birds; i++) {
		ptr = std::make_shared<Bird>(this);
		bird_list.push_back(ptr);
		x[i] = ptr->get_x();
		y[i] = ptr->get_y();
		theta[i] = ptr->get_angle();
	}

	compute_neighbours();
	compute_linklist();
}

VicsekSimulation::VicsekSimulation(const std::list< std::shared_ptr<Bird> >& birds, double velocity, double ang_noise, double length, int seed)
: N_birds(0), vel(velocity), eta(ang_noise), L(length), N_cells(std::ceil(length)), time(0) {
	unsigned i, N_cells_sq = N_cells*N_cells;

	exact_cells = !(std::ceil(length) > length);
	gen = std::mt19937(seed);
	ran_theta = std::uniform_real_distribution<double>(-eta/2, eta/2);
	ran_u = std::uniform_real_distribution<double>(0, 1);

	for (i = 0; i < N_cells_sq; i++) {
		neighbours.push_back(std::vector<unsigned>(8));
		bird_linklist.push_back(std::list< std::shared_ptr<Bird> >());
	}

	for (const std::shared_ptr<Bird>& ptr : birds) {
		bird_list.push_back(ptr);
		ptr->set_sim(this);
		ptr->set_pos(ptr->get_x(), ptr->get_y());
		N_birds++;
	}

	x = Eigen::VectorXd(N_birds);
	y = Eigen::VectorXd(N_birds);
	theta = Eigen::VectorXd(N_birds);

	i = 0;
	for (const std::shared_ptr<Bird>& ptr : birds) {
		x[i] = ptr->get_x();
		y[i] = ptr->get_y();
		theta[i] = ptr->get_angle();
		i++;
	}

	compute_neighbours();
	compute_linklist();
}

void VicsekSimulation::set_N(unsigned N) {
	unsigned i, ind;
	std::shared_ptr<Bird> ptr;
	std::list< std::shared_ptr<Bird> >::iterator it;

	x.conservativeResize(N);
	y.conservativeResize(N);
	theta.conservativeResize(N);

	if (N < N_birds) {
		for (i = N; i < N_birds; i++) {
			ptr = bird_list.back();
			bird_list.pop_back();
			ind = pos_by_indeces(ptr->get_i(), ptr->get_j());
			it = std::find(bird_linklist[ind].begin(), bird_linklist[ind].end(), ptr);
			bird_linklist[ind].erase(it);
		}
	} else {
		for (i = N_birds; i < N; i++) {
			ptr = std::make_shared<Bird>(this);
			bird_list.push_back(ptr);
			bird_linklist[pos_by_indeces(ptr->get_i(), ptr->get_j())].push_back(ptr);
			x[i] = ptr->get_x();
			y[i] = ptr->get_y();
			theta[i] = ptr->get_angle();
		}
	}

	N_birds = N;
}

void VicsekSimulation::set_L(double length) {
	unsigned i, N_cells_new, N_cells_sq_new, N_cells_sq;
	L = length;
	N_cells_new = std::ceil(length);
	exact_cells = !(std::ceil(length) > length);
	bool moreCells = N_cells < N_cells_new;

	if (N_cells == N_cells_new) return;

	N_cells_sq_new = N_cells_new*N_cells_new;
	N_cells_sq = N_cells*N_cells;

	for (i = 0; i < N_cells_sq; i++)
		bird_linklist[i].clear();

	if (moreCells) {
		for (i = N_cells_sq; i < N_cells_sq_new; i++) {
			bird_linklist.push_back(std::list< std::shared_ptr<Bird> >());
			neighbours.push_back(std::vector<unsigned>(8));
		}
	} else {
		bird_linklist.resize(N_cells_sq_new);
		neighbours.resize(N_cells_sq_new);
	}

	N_cells = N_cells_new;

	if (!moreCells) {
		for (const std::shared_ptr<Bird>& ptr : bird_list)
			ptr->set_pos(ptr->get_x(), ptr->get_y());
	}

	compute_linklist();
	compute_neighbours();
}

void VicsekSimulation::push(const std::shared_ptr<Bird>& ptr) {
	N_birds++;

	x.conservativeResize(N_birds);
	y.conservativeResize(N_birds);
	theta.conservativeResize(N_birds);

	bird_list.push_back(ptr);
	ptr->set_sim(this);
	ptr->set_pos(ptr->get_x(), ptr->get_y());
	x[N_birds-1] = ptr->get_x();
	y[N_birds-1] = ptr->get_y();
	theta[N_birds-1] = ptr->get_angle();
	bird_linklist[pos_by_indeces(ptr->get_i(), ptr->get_j())].push_back(ptr);
}

void VicsekSimulation::push(const std::list< std::shared_ptr<Bird> >& birds) {
	unsigned i, n = 0;

	for (const std::shared_ptr<Bird>& ptr : birds) {
		bird_list.push_back(ptr);
		ptr->set_sim(this);
		ptr->set_pos(ptr->get_x(), ptr->get_y());
		bird_linklist[pos_by_indeces(ptr->get_i(), ptr->get_j())].push_back(ptr);
		n++;
	}

	i = N_birds;
	N_birds += n;

	x.conservativeResize(N_birds);
	y.conservativeResize(N_birds);
	theta.conservativeResize(N_birds);

	for (const std::shared_ptr<Bird>& ptr : birds) {
		x[i] = ptr->get_x();
		y[i] = ptr->get_y();
		theta[i] = ptr->get_angle();
	}

}

void VicsekSimulation::randomize() {
	unsigned i;

	for (std::list< std::shared_ptr<Bird> >& list : bird_linklist)
		list.clear();

	i = 0;
	for (const std::shared_ptr<Bird>& ptr : bird_list) {
		ptr->randomize();
		bird_linklist[pos_by_indeces(ptr->get_i(), ptr->get_j())].push_back(ptr);
		x[i] = ptr->get_x();
		y[i] = ptr->get_y();
		theta[i] = ptr->get_angle();
		i++;
	}
}

void VicsekSimulation::change_bird_in_linklist(Bird* b, unsigned i, unsigned j) {
	std::list< std::shared_ptr<Bird> >::iterator it;
	const unsigned ind = pos_by_indeces(i, j);

	if (ind >= N_cells*N_cells) return;

	it = std::find_if(bird_linklist[ind].begin(), bird_linklist[ind].end(), [&](const std::shared_ptr<Bird>& ptr){ return ptr.get() == b; });
	if (it != bird_linklist[ind].end()) {
		// We add the pointer to the new position
		bird_linklist[pos_by_indeces(b->get_i(), b->get_j())].push_back(*it);
		// We delete de pointer from the old position
		bird_linklist[ind].erase(it);
	}
}

void VicsekSimulation::update(unsigned N_steps) {
	unsigned i;

	for (i = 0; i < N_steps; i++) {
		for (const std::shared_ptr<Bird>& bird : bird_list) bird->update_to_new(neighbours, bird_linklist);
		for (const std::shared_ptr<Bird>& bird : bird_list) bird->update_from_new(bird_linklist);
	}
	time += N_steps;
	// Save the information to the vectors
	i = 0;
	for (const std::shared_ptr<Bird>& ptr : bird_list) {
		x[i] = ptr->get_x();
		y[i] = ptr->get_y();
		theta[i] = ptr->get_angle();
		i++;
	}
}

double VicsekSimulation::compute_va() const {
	double sum_sin = .0, sum_cos = .0;
	double ang;

	for (const std::shared_ptr<Bird>& bird : bird_list) {
		ang = bird->get_angle();
		sum_sin += std::sin(ang);
		sum_cos += std::cos(ang);
	}

	return std::sqrt(sum_sin*sum_sin + sum_cos*sum_cos)/N_birds;
}

std::ostream& operator<<(std::ostream& os, const VicsekSimulation& sim) {
	for (const std::shared_ptr<Bird>& bird : sim.bird_list)
		os << sim.time << '\t' << bird->get_x() << '\t' << bird->get_y() << '\t' << bird->get_angle() << '\n';
	return os;
}
