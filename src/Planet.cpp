/*
 * Planet.cpp
 *
 *  Created on: May 29, 2020
 *      Author: rodney
 */

#include "Planet.h"
#include "constants.h"

Planet::Planet(double m, double r)
{
	mass = m;
	radius = r;
	k_g = -(m * G);
}
