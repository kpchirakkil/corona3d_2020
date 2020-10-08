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
	Background_Species(int num_parts, string config_files[], Planet p, double ref_T, double ref_h, string temp_profile_filename, string dens_profile_filename, double profile_bottom, double profile_top);
	virtual ~Background_Species();
	bool check_collision(Particle* p, double dt);
	int get_num_collisions();
	Particle* get_collision_target();
	double get_collision_theta();

private:
	bool use_temp_profile;       // flag for whether or not temperature profile is available
	bool use_dens_profile;       // flag for whether or not density profile is available
	int num_species;             // number of background species in atmosphere
	int num_collisions;          // tracks total number of collisions during simulation
	int collision_target;        // index of particle in bg_parts to be used for next collision
	Planet my_planet;            // contains planet mass, radius, and gravitational constant
	double collision_theta;      // angle (in radians) to be used for next collision
	double ref_temp;             // temperature (in Kelvin) at reference height
	double ref_height;           // height (in cm above planet surface) to extrapolate densities from if no profile available
	double ref_g;                // acceleration due to gravity (G*M/r^2) at reference height
	Distribution_MB* my_dist;    // distribution to be used for initialization of bg particles
	vector<Particle*> bg_parts;  // array of pointers to child particle classes

	double profile_bottom_alt;            //altitude above surface (cm) where atmospheric profiles begin
	double profile_top_alt;               //altitude above surface (cm) where atmospheric profiles end
	vector<double> temp_alt_bins;         // altitude bins from imported temperature profile
	vector<double> Tn;                    // neutral species temperature profile
	vector<double> Ti;                    // ionic species temperature profile
	vector<double> Te;                    // electron temperature profile
	vector<double> dens_alt_bins;         // array of altitude bins imported along with densities
	vector<vector<double>> bg_densities;  // array of densities for each particle in bg_parts
	vector<double> bg_sigma_defaults;     // array of default total cross sections for each particle
	vector<vector<vector<double>>> bg_sigma_tables;  // lookup tables for total cross sections
	vector<double> bg_scaleheights;       // array of scale heights for each particle type
	vector<vector <double>> bg_avg_v;     // array of average (thermal) velocities for each particle

	vector<vector<double>> diff_sigma_energies;
	vector<vector<vector<vector<double>>>> diff_sigma_CDFs;

	// returns collision energy in eV between particle 1 and particle 2
	double calc_collision_e(Particle* p1, Particle* p2);

	// calculates new density of background particle based on radial position and scale height
	double calc_new_density(double ref_density, double scale_height, double r_moved);

	// scans imported differential scattering CDF for new collision theta
	double find_new_theta(int part_index, double energy);

	// get density from imported density profile if available
	double get_density(double alt, vector<double> &dens, double targ_mass);

	// make a new differential cross section CDF and store at diff_sigma_CDFs[index]
	void make_new_CDF(int part_index, int energy_index, vector<double> &angle, vector<double> &sigma);

	//subroutine to set particle type
	Particle* set_particle_type(string type);
};

#endif /* BACKGROUND_SPECIES_HPP_ */
