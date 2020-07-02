/*
 * Background_Species.hpp
 *
 *  Created on: Jun 29, 2020
 *      Author: rodney
 */

#ifndef BACKGROUND_SPECIES_HPP_
#define BACKGROUND_SPECIES_HPP_

#include "Particle_CO.hpp"
#include "Particle_CO2.hpp"
#include "Particle_H.hpp"
#include "Particle_N2.hpp"
#include "Particle_O.hpp"

class Background_Species {
public:
	Background_Species();
	virtual ~Background_Species();
	bool check_collision(double r, double r_ref, double v, double dt);
	int get_num_collisions();
	double get_targ_mass();
	double get_targ_vx();
	double get_targ_vy();
	double get_targ_vz();
	double get_collision_theta();
	void init(double T_bg);

private:
	int num_collisions;
	Particle_CO bg_CO;
	Particle_CO2 bg_CO2;
	Particle_H bg_H;
	Particle_N2 bg_N2;
	Particle_O bg_O;
	double CO_density;
	double CO_sigma;
	double CO_scaleheight;
	double CO_vavg;
	double CO2_density;
	double CO2_sigma;
	double CO2_scaleheight;
	double CO2_vavg;
	double H_density;
	double H_sigma;
	double H_scaleheight;
	double H_vavg;
	double N2_density;
	double N2_sigma;
	double N2_scaleheight;
	double N2_vavg;
	double O_density;
	double O_sigma;
	double O_scaleheight;
	double O_vavg;
	double targ_mass;
	double targ_vx;
	double targ_vy;
	double targ_vz;
	double collision_theta;

	double calc_new_density(double ref_density, double scale_height, double r_moved);

	// returns uniformly distributed random number between 0 and 1
	double get_rand() {return ((double)rand() / RAND_MAX);}
};

#endif /* BACKGROUND_SPECIES_HPP_ */
