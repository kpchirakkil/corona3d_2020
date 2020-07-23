/*
 * Distribution_Hot_O.cpp
 *
 *  Created on: Jul 22, 2020
 *      Author: rodney
 */

#include "Distribution_Hot_O.hpp"

Distribution_Hot_O::Distribution_Hot_O(Planet my_p, double ref_h, double ref_T)
	: Distribution(my_p, ref_h, ref_T) {

	DR160 = 600e6;
	H_DR = 22e3;
	T_ion = 400.0;
	T_e = 1600.0;
	m_ion = 31.9983*constants::amu;

}

Distribution_Hot_O::~Distribution_Hot_O() {

}

void Distribution_Hot_O::init(Particle* p)
{

}
