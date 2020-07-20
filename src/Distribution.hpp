/*
 * Distribution.hpp
 *
 *  Created on: Jul 18, 2020
 *      Author: rodney
 */

#ifndef DISTRIBUTION_HPP_
#define DISTRIBUTION_HPP_

#include "Particle.hpp"

class Distribution {
public:
	Distribution();
	virtual ~Distribution();
	virtual void init(Particle* p) = 0;
};

#endif /* DISTRIBUTION_HPP_ */
