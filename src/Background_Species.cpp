/*
 * Background_Species.cpp
 *
 *  Created on: Jun 29, 2020
 *      Author: rodney
 */

#include "Background_Species.hpp"
#include "constants.hpp"
#include <iostream>

Background_Species::Background_Species() {
	num_collisions = 0;
	CO_density = 0.0;
	CO_sigma = 0.0;
	CO_scaleheight = 0.0;
	CO_vavg = 0.0;
	CO2_density = 0.0;
	CO2_sigma = 0.0;
	CO2_scaleheight = 0.0;
	CO2_vavg = 0.0;
	H_density = 0.0;
	H_sigma = 0.0;
	H_scaleheight = 0.0;
	H_vavg = 0.0;
	N2_density = 0.0;
	N2_sigma = 0.0;
	N2_scaleheight = 0.0;
	N2_vavg = 0.0;
	O_density = 0.0;
	O_sigma = 0.0;
	O_scaleheight = 0.0;
	O_vavg = 0.0;
	targ_mass = 0.0;
	targ_vx = 0.0;
	targ_vy = 0.0;
	targ_vz = 0.0;
	collision_theta = 0.0;
}

Background_Species::~Background_Species() {

}

double Background_Species::calc_new_density(double ref_density, double scale_height, double r_moved)
{
	return ref_density*exp(r_moved/scale_height);
}

// check to see if a collision occurred and set targ values if so
bool Background_Species::check_collision(double r, double r_ref, double v, double dt)
{
	// get densities at current location
	double r_moved = r_ref - r;
	double CO_temp_density = calc_new_density(CO_density, CO_scaleheight, r_moved);
	double CO2_temp_density = calc_new_density(CO2_density, CO2_scaleheight, r_moved);
	double H_temp_density = calc_new_density(H_density, H_scaleheight, r_moved);
	double N2_temp_density = calc_new_density(N2_density, N2_scaleheight, r_moved);
	double O_temp_density = calc_new_density(O_density, O_scaleheight, r_moved);

	// determine if test particle collided
	double u = get_rand();
	double tau = v*dt*(O_sigma*O_temp_density +
			           CO_sigma*CO_temp_density +
					   N2_sigma*N2_temp_density +
					   CO2_sigma*CO2_temp_density +
					   H_sigma*H_temp_density);

	if (u > exp(-tau))
	{
		num_collisions++;

		// pick target species mass, set targ values
		u = get_rand();
		if (u < O_temp_density/(O_temp_density+CO_temp_density+N2_temp_density+CO2_temp_density))
		{
			targ_mass = bg_O.get_mass();  // collision with O
			bg_O.init_particle_vonly_MB(O_vavg);
			targ_vx = bg_O.get_vx();
			targ_vy = bg_O.get_vy();
			targ_vz = bg_O.get_vz();
		}
		else if (u < (O_temp_density+N2_temp_density)/(O_temp_density+CO_temp_density+N2_temp_density+CO2_temp_density))
		{
			targ_mass = bg_N2.get_mass();  //collision with N2
			bg_N2.init_particle_vonly_MB(N2_vavg);
			targ_vx = bg_N2.get_vx();
			targ_vy = bg_N2.get_vy();
			targ_vz = bg_N2.get_vz();
		}
		else if (u < (O_temp_density+CO_temp_density+N2_temp_density)/(O_temp_density+CO_temp_density+N2_temp_density+CO2_temp_density))
		{
			targ_mass = bg_CO.get_mass();     // collision with CO
			bg_CO.init_particle_vonly_MB(CO_vavg);
			targ_vx = bg_CO.get_vx();
			targ_vy = bg_CO.get_vy();
			targ_vz = bg_CO.get_vz();
		}
		else
		{
			targ_mass = bg_CO2.get_mass();    // collision with CO2
			bg_CO2.init_particle_vonly_MB(CO2_vavg);
			targ_vx = bg_CO2.get_vx();
			targ_vy = bg_CO2.get_vy();
			targ_vz = bg_CO2.get_vz();
		}

		return true;
	}
	else
	{
		return false;
	}
}

int Background_Species::get_num_collisions()
{
	return num_collisions;
}

double Background_Species::get_targ_mass()
{
	return targ_mass;
}

double Background_Species::get_targ_vx()
{
	return targ_vx;
}

double Background_Species::get_targ_vy()
{
	return targ_vy;
}

double Background_Species::get_targ_vz()
{
	return targ_vz;
}

double Background_Species::get_collision_theta()
{
	return collision_theta;
}

void Background_Species::init(double T_bg)
{
	num_collisions = 0;
	CO_density = 1.0e13;
	CO_sigma = 1.85e-18;
	CO_scaleheight = constants::k_b*T_bg/(bg_CO.get_mass()*3.31);
	CO_vavg = sqrt(constants::k_b*T_bg/bg_CO.get_mass());
	CO2_density = 6.68e13;
	CO2_sigma = 2.0e-18;
	CO2_scaleheight = constants::k_b*T_bg/(bg_CO2.get_mass()*3.31);
	CO2_vavg = sqrt(constants::k_b*T_bg/bg_CO2.get_mass());
	H_density = 0.0;
	H_sigma = 0.0;
	H_scaleheight = constants::k_b*T_bg/(bg_H.get_mass()*3.31);
	H_vavg = sqrt(constants::k_b*T_bg/bg_H.get_mass());
	N2_density = 3.1e13;
	N2_sigma = 1.85e-18;
	N2_scaleheight = constants::k_b*T_bg/(bg_N2.get_mass()*3.31);
	N2_vavg = sqrt(constants::k_b*T_bg/bg_N2.get_mass());
	O_density = 2.64e13;
	O_sigma = 6.4e-19;
	O_scaleheight = constants::k_b*T_bg/(bg_O.get_mass()*3.31);
	O_vavg = sqrt(constants::k_b*T_bg/bg_O.get_mass());
}
