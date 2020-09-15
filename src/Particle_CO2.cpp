/*
 * Particle_CO2.cpp
 *
 *  Created on: Jun 23, 2020
 *      Author: rodney
 */

#include "Particle_CO2.hpp"

const double Particle_CO2::mass = 44.0095*constants::amu;

Particle_CO2::Particle_CO2() {

}

Particle_CO2::~Particle_CO2() {

}

double Particle_CO2::get_mass()
{
	return mass;
}
