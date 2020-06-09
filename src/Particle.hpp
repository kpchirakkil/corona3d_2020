/*
 * Particle.h
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */

#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_

#include <cstdlib>

class Particle {
public:
	Particle();
	void initParticle_MB(double r, double v_avg); // init particle from MB distribution

	bool active;            // flag for whether particle is active, i.e. should still be considered in the simulation
	double radius;          // radius from center of planet	[m]
	double inverseRadius;   // inverse radius (for computational efficiency) [m^-1]
	double position[3];     // position vector [m,m,m]
	double velocity[3];     // velocity vector [m/s,m/s,m/s]

private:
	// returns uniformly distributed random number between 0 and 1
	double getRand() {return ((double)rand() / RAND_MAX);}
};

#endif /* PARTICLE_HPP_ */
