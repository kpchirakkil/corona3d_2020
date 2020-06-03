/*
 * Planet.cpp
 *
 *  Created on: May 29, 2020
 *      Author: rodney
 */

#include "Planet.h"
#include "constants.h"

Planet::Planet(double m, double r, double T, double part_m)
{
	mass = m;
	radius = r;
	k_g = -(m * G);
	T_bg = T;
	particle_m = part_m;
}
