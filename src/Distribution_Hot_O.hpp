/*
 * Distribution_Hot_O.hpp
 *
 *  Created on: Jul 22, 2020
 *      Author: rodney
 */

#ifndef DISTRIBUTION_HOT_O_HPP_
#define DISTRIBUTION_HOT_O_HPP_

#include "Distribution.hpp"

class Distribution_Hot_O: public Distribution {
public:
	Distribution_Hot_O(Planet my_p, double ref_h, double ref_T);
	virtual ~Distribution_Hot_O();
	void init(Particle* p);

private:
	double DR160;  // [m^3/s] dissociative recombination rate at 200km
	double H_DR;   // [m] dissociative recombination scale height
	double T_ion;  // [K] ion temperature
	double T_e;    // [K] electron temperature
	double m_ion;  // [kg] O2+ ion mass

	double E_rot(double B, double T);
	void gen_mb(double vavg, double v[]);
};

#endif /* DISTRIBUTION_HOT_O_HPP_ */
