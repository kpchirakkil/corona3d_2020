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
	double m_H;                // [g] neutral H mass
	double m_Hion;             // [g] H+ ion mass
	double H_Hplus_rate_coeff; // [cm^3/s] rate coefficient for H+ + H -> H* + H+
	vector<vector<double>> H_profile;
	vector<vector<double>> Hplus_profile;
	vector<vector<double>> temp_profile;
	vector<vector<double>> H_Hplus_CDF;

	// scans H_Hplus_CDF for new particle radius
	double get_new_radius_H_Hplus();

	// generate H_Hplus_CDF for given altitude range using imported density/temp profiles
	void make_H_Hplus_CDF(double lower_alt, double upper_alt);
};

#endif /* DISTRIBUTION_HOT_H_HPP_ */
