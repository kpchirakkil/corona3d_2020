/*
 * Particle.cpp
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */

#include "Particle.hpp"
#include "constants.hpp"

Particle::Particle()
{
	active = true;
	radius = 0.0;
	inverseRadius = 0.0;
	position[0] = position[1] = position[2] = 0.0;
	velocity[0] = velocity[1] = velocity[2] = 0.0;
}

// initialize a single particle at given radius using Maxwell-Boltzmann avg v
void Particle::initParticle_MB(double r, double v_avg)
{
	radius = r;
	double phi = twopi*(getRand());
	double u = 2.0*getRand() - 1;
	inverseRadius = 1.0/r;
	position[0] = r*sqrt(1-(u*u))*cos(phi);
	position[1] = r*sqrt(1-(u*u))*sin(phi);
	position[2] = r*u;

	double randnum1 = getRand();
	double randnum2 = getRand();
	double randnum3 = getRand();
	double randnum4 = getRand();

	randnum1 = v_avg*sqrt(-2.0*log(1.0-randnum1));
	randnum2 = twopi*randnum2;
	randnum3 = v_avg*sqrt(-2.0*log(1.0-randnum3));
	randnum4 = twopi*randnum4;

	velocity[0] = randnum1*cos(randnum2);
	velocity[1] = randnum1*sin(randnum2);
	velocity[2] = randnum3*cos(randnum4);
}
