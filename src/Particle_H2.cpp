/*
 * Particle_H2.cpp
 *
 *  Created on: July 16, 2024
 *      Author: Grace
 */

#include "Particle_H2.hpp"

const double Particle_H2::mass = 2.016*constants::amu;
const string Particle_H2::name = "H2";

Particle_H2::Particle_H2() {

}

Particle_H2::~Particle_H2() {

}

double Particle_H2::get_mass() const
{
	return mass;
}

string Particle_H2::get_name() const
{
	return name;
}
