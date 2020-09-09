/*
 * Atmosphere.hpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#ifndef ATMOSPHERE_HPP_
#define ATMOSPHERE_HPP_

#include <vector>
#include <iomanip>
#include "Background_Species.hpp"
#include "Distribution_Hot_H.hpp"
#include "Distribution_Hot_O.hpp"
#include "Distribution_Import.hpp"
#include "Distribution_MB.hpp"
#include "Common_Functions.hpp"
using namespace std;

class Atmosphere {
public:
	Atmosphere(int n, Planet p, vector<Particle*> parts, Distribution* dist, Background_Species bg, double T, double ref_h, string temp_profile);
	virtual ~Atmosphere();

	void output_positions(std::string datapath);
	void output_altitude_distro(double bin_width, int num_bins, std::string datapath);
	void output_velocity_distro(double bin_width, int num_bins, std::string datapath);
	void run_simulation(double dt, int num_steps);

private:
	Common_Functions common;
	int num_parts;                      // number of particles initially spawned
	int active_parts;                   // number of active particles
	Planet my_planet;                   // contains planet mass and radius
	vector<Particle*> my_parts;         // array of particles to be tracked
	Distribution* my_dist;              // distribution class to initialize particles
	Background_Species bg_species;      // background species used for collisions
	double T_bg;                        // [K] background temp where simulation starts
	double ref_height;                  // [cm] altitude above planet surface of model bottom

	vector<double> alt_bins;
	vector<double> Tn;
	vector<double> Ti;
	vector<double> Te;
};

#endif /* ATMOSPHERE_HPP_ */
