/*
 * Atmosphere.hpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#ifndef ATMOSPHERE_HPP_
#define ATMOSPHERE_HPP_

#include <vector>
#include "Planet.hpp"
#include "Particle_H.hpp"

class Atmosphere {
public:
	Atmosphere(); // initialize with defaults (10K MB-distributed H atoms 160km above Venus)
	Atmosphere(int n, Planet p, double T, double model_b);
	virtual ~Atmosphere();

	void do_timestep(double dt);
	void output_positions(std::string datapath);
	void output_velocity_distro(double bin_width, int num_bins, std::string datapath);
	void run_simulation(double dt, int numSteps);

private:
	int N;                        // number of particles to track
	Planet myPlanet;              // contains planet mass and radius
	std::vector<Particle_H> myH;  // array of H atoms to be tracked
	double T_bg;                  // [K] background temp where simulation starts
	double model_bottom;          // [m] altitude above planet surface of model bottom
};

#endif /* ATMOSPHERE_HPP_ */
