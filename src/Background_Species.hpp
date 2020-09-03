/*
 * Background_Species.hpp
 *
 *  Created on: Jun 29, 2020
 *      Author: rodney
 */

#ifndef BACKGROUND_SPECIES_HPP_
#define BACKGROUND_SPECIES_HPP_

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include "Particle_CO.hpp"
#include "Particle_CO2.hpp"
#include "Particle_H.hpp"
#include "Particle_N2.hpp"
#include "Particle_O.hpp"
#include "Distribution_MB.hpp"
#include "Planet.hpp"
#include "Common_Functions.hpp"
using namespace std;

class Background_Species {
public:
	Background_Species();
	Background_Species(int n, Planet p, double T, double h, Distribution_MB* dist, Particle* bg_p[], double bg_d[], double bg_s[], string dens_profile_filename);
	virtual ~Background_Species();
	bool check_collision(double r, double v, double dt, double T);
	int get_num_collisions();
	Particle* get_collision_target();
	double get_collision_theta();
	void import_CDF(string filename);

private:
	Common_Functions common;
	bool use_dens_profile;       // flag for whether or not density profile is available
	int num_species;             // number of background species in atmosphere
	int num_collisions;          // tracks total number of collisions during simulation
	int collision_target;        // index of particle in bg_parts to be used for next collision
	Planet my_planet;            // contains planet mass, radius, and gravitational constant
	double collision_theta;      // angle (in radians) to be used for next collision
	double ref_temp;             // temperature (in Kelvin) at reference height
	double ref_height;           // height (in meters above planet surface) where simulation starts
	double ref_g;                // acceleration due to gravity (G*M/r^2) at reference height
	Distribution_MB* my_dist;    // distribution to be used for initialization of bg particles
	vector<Particle*> bg_parts;  // array of pointers to child particle classes
	vector<double> alt_bins;     // array of altitude bins imported along with densities
	vector<vector<double>> bg_densities;  // array of densities for each particle in bg_parts
	vector<double> bg_sigmas;             // array of total cross sections for each particle
	vector<double> bg_scaleheights;       // array of scale heights for each particle type
	vector<double> bg_avg_v;              // array of average (thermal) velocities for each particle
	Matrix<double, 180, 2> cdf;           // 180x2 array of cos(theta) vs diff. cross section CDF values

	// calculates new density of background particle based on radial position and scale height
	double calc_new_density(double ref_density, double scale_height, double r_moved);

	// scans imported differential scattering CDF for new collision theta
	double find_new_theta();

	// get density from imported density profile if available
	double get_density(double r, vector<double> &dens);
};

#endif /* BACKGROUND_SPECIES_HPP_ */
