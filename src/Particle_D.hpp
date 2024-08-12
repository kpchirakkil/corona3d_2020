/*
 * Particle_D.hpp
 *
 *  Created on: Aug 2, 2024
 *      Author: Grace
 */

#ifndef PARTICLE_D_HPP_
#define PARTICLE_D_HPP_

#include "Particle.hpp"

class Particle_D: public Particle {
public:
	Particle_D();
	virtual ~Particle_D();
	static const double mass;
	static const string name;
	double get_mass() const;
	string get_name() const;
};

#endif /* PARTICLE_D_HPP_ */
