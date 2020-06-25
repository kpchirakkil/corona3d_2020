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
#include "Particle_H.hpp"
using namespace std;

int main(int argc, char* argv[])
{
	int timesteps = 1000;
	double dt = 0.5;
	Planet venus;
	venus.init();
	Atmosphere myAtmosphere(20000, venus, 350.0, 160e3);

	myAtmosphere.output_velocity_distro(100.0, 150, "/home/rodney/Documents/coronaTest/vdist.out");

	for (int i=0; i<timesteps; i++)
	{
		myAtmosphere.output_positions("/home/rodney/Documents/coronaTest/data/positions" + to_string(i) + ".out");
		myAtmosphere.do_timestep(dt);
	}

	myAtmosphere.output_velocity_distro(100.0, 150, "/home/rodney/Documents/coronaTest/vdist2.out");

	return 0;
}
