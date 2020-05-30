/*
 * phys.h
 *
 *  Created on: May 29, 2020
 *      Author: rodney
 */

#ifndef PHYS_H_
#define PHYS_H_

#include "Particle.h"
#include "Planet.h"

void initPlanet(Planet plnt);
void initParticles(Planet plnt, Particle parts[], int n);
void stepParticles(Particle parts[], double dt, int n, double k_g);
double getRand();

#endif /* PHYS_H_ */
