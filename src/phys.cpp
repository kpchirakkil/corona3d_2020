/*
 * phys.cpp
 *
 *  Created on: May 29, 2020
 *      Author: rodney
 */

#include "phys.h"
#include "constants.h"
#include "Particle.h"
#include "Planet.h"
#include <cstdlib>
#include <cmath>

void initPlanet(Planet plnt)
{

}

void initParticles(Planet plnt, Particle parts[], int n)
{
	double scaleHeight = 20e3; // should calc/assign in planet class eventually
	for (int i=0; i<n; i++)
	{
		double r = plnt.radius + 160e3 - log(getRand())*scaleHeight;
		double phi = twopi*(getRand());
		double u = 2.0*getRand() - 1;
		parts[i].radius = r;
		parts[i].inverseRadius = 1.0/r;
		parts[i].position[0] = r*sqrt(1-(u*u))*cos(phi);
		parts[i].position[1] = r*sqrt(1-(u*u))*sin(phi);
		parts[i].position[2] = r*u;
	}
}

double getRand()
{
	return (double) rand() / RAND_MAX;
}
