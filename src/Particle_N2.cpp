/*
 * Particle_N2.cpp
 *
 *  Created on: Jun 23, 2020
 *      Author: rodney
 */

#include "Particle_N2.hpp"
#include "constants.hpp"

const double Particle_N2::mass = 28.0134*constants::amu;

Particle_N2::Particle_N2() {

}

Particle_N2::~Particle_N2() {

}

double Particle_N2::get_mass()
{
	return mass;
}
