/*
 * Atmosphere.hpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#ifndef ATMOSPHERE_HPP_
#define ATMOSPHERE_HPP_

#include <vector>
#include "Background_Species.hpp"

class Atmosphere {
public:
	Atmosphere(int n, Planet p, Background_Species bg, double T, double ref_h);
	virtual ~Atmosphere();

	void output_positions(std::string datapath);
	void output_velocity_distro(double bin_width, int num_bins, std::string datapath);
	void run_simulation(double dt, int num_steps);

private:
	int N;                              // number of particles to track
	Planet my_planet;                   // contains planet mass and radius
	std::vector<Particle_O> my_parts;   // array of atoms to be tracked
	Background_Species bg_species;      // background species used for collisions
	double T_bg;                        // [K] background temp where simulation starts
	double ref_height;                  // [m] altitude above planet surface of model bottom
};

#endif /* ATMOSPHERE_HPP_ */
