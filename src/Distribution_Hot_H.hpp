/*
 * Distribution_Hot_H.hpp
 *
 *  Created on: Jul 31, 2020
 *      Author: rodney
 */

#ifndef DISTRIBUTION_HOT_H_HPP_
#define DISTRIBUTION_HOT_H_HPP_

#include "Distribution.hpp"

class Distribution_Hot_H: public Distribution {
public:
	Distribution_Hot_H(Planet my_p, double ref_h, double ref_T);
	virtual ~Distribution_Hot_H();
	void init(Particle* p);

private:
	double T_ion;              // [K] ion temperature
	double T_e;                // [K] electron temperature
	double m_Hion;             // [g] H+ ion mass
	double H_Hplus_rate_coeff; // [cm^3/s] rate coefficient for H+ + H -> H* + H+

	vector<vector<double>> H_profile;
	vector<vector<double>> Hplus_profile;
	vector<vector<double>> temp_profile;
	vector<double> alt_bins;
	vector<double> H_Hplus_CDF;

	// scans H_Hplus_CDF for new particle radius
	double get_new_radius_H_Hplus();
};

#endif /* DISTRIBUTION_HOT_H_HPP_ */
