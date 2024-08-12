/*
 * Particle_D.cpp
 *
 *  Created on: Aug 2, 2024
 *      Author: Grace
 */

#include "Particle_D.hpp"

const double Particle_D::mass = 2.014*constants::amu;
const string Particle_D::name = "D";

Particle_D::Particle_D() {

}

Particle_D::~Particle_D() {
}

double Particle_D::get_mass() const
{
	return mass;
}

string Particle_D::get_name() const
{
	return name;
}
