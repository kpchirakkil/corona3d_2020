/*
 * Planet.cpp
 *
 *  Created on: May 29, 2020
 *      Author: rodney
 */

#include "Planet.hpp"
#include "constants.hpp"

Planet::Planet()
{
	mass = 0.0;
	radius = 0.0;
	k_g = 0.0;
}

void Planet::init()
{
	mass = 4.8675e24;
	radius = 6051.8e3;
	k_g = -(mass * G);
}

void Planet::init(double m, double r)
{
	mass = m;
	radius = r;
	k_g = -(mass * G);
}
