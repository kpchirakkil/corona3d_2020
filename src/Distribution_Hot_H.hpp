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
	void init(shared_ptr<Particle> p);

private:
	double m_H;                // [g] neutral H mass
	double m_Hplus;            // [g] H+ ion mass
	double m_HCOplus;          // [g] HCO+ ion mass
	double m_CO;               // [g] CO mass
	double H_Hplus_rate_coeff; // [cm^3/s] rate coefficient for H+ + H -> H* + H+
	double HCOplus_DR_rate_coeff;  // [cm^3/s] rate coefficient for HCO+ + e -> H* + CO
	vector<vector<double>> H_profile;
	vector<vector<double>> Hplus_profile;
	vector<vector<double>> HCOplus_profile;
	vector<vector<double>> electron_profile;
	vector<vector<double>> temp_profile;
	vector<vector<double>> H_Hplus_CDF;
	vector<vector<double>> HCOplus_DR_CDF;

	// scans HCOplus_DR_CDF for new particle radius
	double get_new_radius_HCOplus_DR();

	// scans H_Hplus_CDF for new particle radius
	double get_new_radius_H_Hplus();

	// init particle using H_Hplus mechanism
	void init_H_Hplus_particle(shared_ptr<Particle> p);

	// init particle using HCOplus_DR mechanism
	void init_HCOplus_DR_particle(shared_ptr<Particle> p);

	// generate HCOplus_DR_CDF for given altitude range using imported density/temp profiles
	void make_HCOplus_DR_CDF(double lower_alt, double upper_alt);

	// generate H_Hplus_CDF for given altitude range using imported density/temp profiles
	void make_H_Hplus_CDF(double lower_alt, double upper_alt);
};

#endif /* DISTRIBUTION_HOT_H_HPP_ */
