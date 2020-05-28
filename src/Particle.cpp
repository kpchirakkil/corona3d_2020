/*
 * Particle.cpp
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */

#include "Particle.h"

Particle::Particle() {
	active = true;
	radius = 0.0;
	inverseRadius = 0.0;
	position[0] = position[1] = position[2] = 0.0;
	velocity[0] = velocity[1] = velocity[2] = 0.0;
}

Particle::~Particle() {
	// TODO Auto-generated destructor stub
}

