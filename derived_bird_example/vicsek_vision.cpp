#include "vicsek_vision.hpp"

void BirdVision::sum_birds_in_view(const std::list< std::shared_ptr<Bird> >& birds, double& sum_sin, double& sum_cos) const {
	double dif_ang;
	double bird_angle = get_angle();

	for (const std::shared_ptr<Bird>&  p : birds) {
		if (p.get() == this) {
			// If the bird is "interacting" with itself we just sum the contribution
			sum_sin += std::sin(bird_angle);
			sum_cos += std::cos(bird_angle);
		} else {
			// Gets the angular difference between the direction of the bird
			// and the angle with the other bird
			dif_ang = get_angle_to(p) - bird_angle;
			// If this difference is bigger than pi or smaller than -pi
			// we can do better (the maximum distance between to directions is pi)
			while (dif_ang > std::numbers::pi) dif_ang -= 2*std::numbers::pi;
			while (dif_ang < -std::numbers::pi) dif_ang += 2*std::numbers::pi;

			// We check if we have to sum the contribution of the bird
			if (std::abs(dif_ang) < ang_field_view*0.5 && distance_sq(p) < 1) {
				sum_sin += std::sin(p->get_angle());
				sum_cos += std::cos(p->get_angle());
			}
		}
	}
}

double BirdVision::get_angle_to(const std::shared_ptr<Bird>& pBird) const {
	double L = get_sim_ptr()->get_L();
	double dx = pBird->get_x() - get_x();
	double dy = pBird->get_y() - get_y();


	// We check the periodic boundary conditions
	if (dx > 0.5*L) dx -= L;
	else if (dx < -0.5*L) dx += L;

	if (dy > 0.5*L) dy -= L;
	else if (dy < -0.5*L) dy += L;

	return std::atan2(dy, dx);
}

double BirdVision::average_angle(const std::vector< std::vector<unsigned> >& neighbours, const std::vector< std::list< std::shared_ptr<Bird> > >& link_list) const {
	unsigned i = get_i(), j = get_j();
	VicsekSimulation* sim = get_sim_ptr();

	unsigned ind = sim->pos_by_indeces(i, j);
	double sum_sin = .0, sum_cos = .0;
	unsigned n, N_neigh = neighbours[0].size(), N_cells;

	sum_birds_in_view(link_list[ind], sum_sin, sum_cos);

	for (n = 0; n < N_neigh; n++) {
		sum_birds_in_view(link_list[neighbours[ind][n]], sum_sin, sum_cos);
	}

	// Some additional checks if the number of cells does not exactly match the length
	if (!sim->get_exact_cells()) {
		N_cells = sim->get_N_cells();
		if (i == 0) {
			n = sim->pos_by_indeces(N_cells-1, j);
			sum_birds_in_view(link_list[neighbours[n][1]], sum_sin, sum_cos);
			sum_birds_in_view(link_list[neighbours[n][4]], sum_sin, sum_cos);
			sum_birds_in_view(link_list[neighbours[n][5]], sum_sin, sum_cos);
		} else if (i == N_cells - 2) {
				n = sim->pos_by_indeces(N_cells-1, j);
				sum_birds_in_view(link_list[neighbours[n][3]], sum_sin, sum_cos);
				sum_birds_in_view(link_list[neighbours[n][6]], sum_sin, sum_cos);
				sum_birds_in_view(link_list[neighbours[n][7]], sum_sin, sum_cos);

		}
		if (j == 0) {
			n = sim->pos_by_indeces(i, N_cells-1);
			sum_birds_in_view(link_list[neighbours[n][0]], sum_sin, sum_cos);
			sum_birds_in_view(link_list[neighbours[n][4]], sum_sin, sum_cos);
			sum_birds_in_view(link_list[neighbours[n][7]], sum_sin, sum_cos);
		} else if (j == N_cells - 2) {
				n = sim->pos_by_indeces(i, N_cells-1);
				sum_birds_in_view(link_list[neighbours[n][2]], sum_sin, sum_cos);
				sum_birds_in_view(link_list[neighbours[n][5]], sum_sin, sum_cos);
				sum_birds_in_view(link_list[neighbours[n][6]], sum_sin, sum_cos);

		}

	}

	// Technically what we should do is the average of the sums
	// but since std::atan2 only cares about the proportion
	// there is no point in doing the averages
	return std::atan2(sum_sin, sum_cos);
}
