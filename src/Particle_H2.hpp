/*
 * Particle_H2.hpp
 *
 *  Created on: July 16, 2024
 *      Author: Grace
 */

#ifndef PARTICLE_H2_HPP_
#define PARTICLE_H2_HPP_

#include "Particle.hpp"

class Particle_H2: public Particle {
public:
	Particle_H2();
	virtual ~Particle_H2();
	static const double mass;
	static const string name;
	double get_mass() const;
	string get_name() const;
};

#endif /* PARTICLE_H2_HPP_ */
