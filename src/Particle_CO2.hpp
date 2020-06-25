/*
 * Particle_CO2.hpp
 *
 *  Created on: Jun 23, 2020
 *      Author: rodney
 */

#ifndef PARTICLE_CO2_HPP_
#define PARTICLE_CO2_HPP_

#include "Particle.hpp"

class Particle_CO2: public Particle {
public:
	Particle_CO2();
	virtual ~Particle_CO2();
	static const double mass;
	double get_mass();
};

#endif /* PARTICLE_CO2_HPP_ */
