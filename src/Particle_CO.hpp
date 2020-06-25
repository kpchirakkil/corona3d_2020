/*
 * Particle_CO.hpp
 *
 *  Created on: Jun 23, 2020
 *      Author: rodney
 */

#ifndef PARTICLE_CO_HPP_
#define PARTICLE_CO_HPP_

#include "Particle.hpp"

class Particle_CO: public Particle {
public:
	Particle_CO();
	virtual ~Particle_CO();
	static const double mass;
	double get_mass();
};

#endif /* PARTICLE_CO_HPP_ */
