/*
 * Distribution.hpp
 *
 *  Created on: Jul 18, 2020
 *      Author: rodney
 */

#ifndef DISTRIBUTION_HPP_
#define DISTRIBUTION_HPP_

#include "Particle.hpp"
#include "Planet.hpp"

class Distribution {
public:
	Distribution(Planet my_p, double ref_h, double ref_T);
	virtual ~Distribution();
	virtual void init(Particle* p) = 0;

protected:
	Planet my_planet;
	double ref_height;
	double ref_radius;
	double ref_temp;

	// sets v_in[] to a Maxwell-Boltzmann velocity vector based on vavg (average velocity)
	void gen_mb(double vavg, double v_in[]);
};

#endif /* DISTRIBUTION_HPP_ */
