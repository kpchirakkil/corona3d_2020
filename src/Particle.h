/*
 * Particle.h
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_

class Particle {
public:
	Particle();
	virtual ~Particle();
//private:
	bool active;			// flag for whether particle is active, i.e. should still be considered in the simulation
	double radius;			// radius from center of planet	[m]
	double inverseRadius;	// inverse radius (for computational efficiency) [m^-1]
	double position[3];		// position vector [m,m,m]
	double velocity[3];		// velocity vector [m/s,m/s,m/s]
};

#endif /* PARTICLE_H_ */
