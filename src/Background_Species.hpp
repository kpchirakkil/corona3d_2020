/*
 * Background_Species.hpp
 *
 *  Created on: Jun 29, 2020
 *      Author: rodney
 */

#ifndef BACKGROUND_SPECIES_HPP_
#define BACKGROUND_SPECIES_HPP_

#include <vector>
#include "Particle_CO.hpp"
#include "Particle_CO2.hpp"
#include "Particle_H.hpp"
#include "Particle_N2.hpp"
#include "Particle_O.hpp"
#include "Planet.hpp"

class Background_Species {
public:
	Background_Species();
	Background_Species(int n, Planet p, double T, double h, Particle* bg_p[], double bg_d[], double bg_s[]);
	virtual ~Background_Species();
	bool check_collision(double r, double v, double dt);
	int get_num_collisions();
	Particle* get_collision_target();
	double get_collision_theta();

private:
	int num_species;
	int num_collisions;
	int collision_target;
	Planet my_planet;
	double collision_theta;
	double ref_temp;
	double ref_height;
	double ref_g;
	std::vector<Particle*> bg_parts;
	std::vector<double> bg_densities;
	std::vector<double> bg_sigmas;
	std::vector<double> bg_scaleheights;
	std::vector<double> bg_avg_v;

	double calc_new_density(double ref_density, double scale_height, double r_moved);

	// returns uniformly distributed random number between 0 and 1
	double get_rand() {return ((double)rand() / (double)RAND_MAX);}
};

#endif /* BACKGROUND_SPECIES_HPP_ */
