/*
 * Particle_O.hpp
 *
 *  Created on: Jun 24, 2020
 *      Author: rodney
 */

#ifndef PARTICLE_O_HPP_
#define PARTICLE_O_HPP_

#include "Particle.hpp"

class Particle_O: public Particle {
public:
	Particle_O();
	virtual ~Particle_O();
	static const double mass;
	double get_mass();
};

#endif /* PARTICLE_O_HPP_ */
