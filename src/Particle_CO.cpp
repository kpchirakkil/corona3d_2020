/*
 * Particle_CO.cpp
 *
 *  Created on: Jun 23, 2020
 *      Author: rodney
 */

#include "Particle_CO.hpp"

const double Particle_CO::mass = 28.0101*constants::amu;

Particle_CO::Particle_CO() {

}

Particle_CO::~Particle_CO() {

}

double Particle_CO::get_mass()
{
	return mass;
}
