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

Planet::~Planet()
{

}

void Planet::init()
{
	mass = 4.8675e24;              // Venus mass by default
	radius = 6051.8e3;             // Venus radius by default
	k_g = -(mass * constants::G);  // planetary gravitational constant
}

void Planet::init(double m, double r)
{
	mass = m;
	radius = r;
	k_g = -(mass * constants::G);
}

double Planet::get_mass()
{
	return mass;
}

double Planet::get_radius()
{
	return radius;
}

double Planet::get_k_g()
{
	return k_g;
}
