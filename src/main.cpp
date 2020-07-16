/*
 * main.cpp
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "constants.hpp"
#include "Planet.hpp"
#include "Atmosphere.hpp"
using namespace std;

int main(int argc, char* argv[])
{
	// simulation parameters
	int num_testparts = 10000;
	int timesteps = 1000;
	double dt = 0.05;
	double ref_alt = 200e3;
	double bg_temp = 277.6;
	double planet_mass = 6.4185e23;
	double planet_radius = 3.397e6;
	Planet mars;
	mars.init(planet_mass, planet_radius);

	// initialize background species
	int num_bgparts = 4;
	Particle* bg_parts[] = {new Particle_O(),
							new Particle_N2(),
							new Particle_CO(),
							new Particle_CO2()};
	double bg_dens[] = {2.64e13, 1.0e13, 3.1e13, 6.68e13};
	double bg_sigs[] = {6.4e-19, 1.85e-18, 1.85e-18, 2.0e-18};
	Background_Species bg_spec(num_bgparts, mars, bg_temp, ref_alt, bg_parts, bg_dens, bg_sigs);

	// initialize atmosphere and run simulation
	Atmosphere my_atmosphere(num_testparts, mars, bg_spec, bg_temp, ref_alt);
	my_atmosphere.output_velocity_distro(100.0, 150, "/home/rodney/Documents/coronaTest/vdist.out");
	my_atmosphere.run_simulation(dt, timesteps);
	my_atmosphere.output_velocity_distro(100.0, 150, "/home/rodney/Documents/coronaTest/vdist2.out");

	return 0;
}
