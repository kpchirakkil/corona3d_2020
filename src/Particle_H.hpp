/*
 * Particle_H.hpp
 *
 *  Created on: Jun 19, 2020
 *      Author: rodney
 */

#ifndef PARTICLE_H_HPP_
#define PARTICLE_H_HPP_

#include "Particle.hpp"

class Particle_H: public Particle {
public:
	Particle_H();
	virtual ~Particle_H();
	static const double mass;
	double get_mass();
};

#endif /* PARTICLE_H_HPP_ */
