/*
 * Particle_O.cpp
 *
 *  Created on: Jun 24, 2020
 *      Author: rodney
 */

#include "Particle_O.hpp"
#include "constants.hpp"

const double Particle_O::mass = 15.9994*constants::amu;

Particle_O::Particle_O() {

}

Particle_O::~Particle_O() {

}

double Particle_O::get_mass()
{
	return mass;
}
