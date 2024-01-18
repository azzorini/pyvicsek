#pragma once

#include <cmath>
#include <vector>
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

	void set_width_view(double ViewWidth) { ang_field_view = ViewWidth; }

	double get_width_view() const { return ang_field_view; }

	// Gets the angle between this bird and the bird given by the pointer
	double get_angle_to(const std::shared_ptr<Bird>&) const;

	// Override the average angle function
	virtual double average_angle(const std::vector< std::vector<unsigned> >&, const std::vector< std::list< std::shared_ptr<Bird> > >&) const;
};
