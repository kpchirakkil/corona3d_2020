/*
 * Distribution_Hot_H.cpp
 *
 *  Created on: Jul 31, 2020
 *      Author: rodney
 */

#include "Distribution_Hot_H.hpp"

Distribution_Hot_H::Distribution_Hot_H(Planet my_p, double ref_h, double ref_T)
	: Distribution(my_p, ref_h, ref_T) {

	m_Hion = 1.00728*constants::amu;
	H_Hplus_rate_coeff = 8.7e-10;

	temp_profile.resize(4);
	H_profile.resize(2);
	Hplus_profile.resize(2);

	string temp_prof_filename = "/home/rodney/git/corona3d_2020/src/inputs/Mars/MarsTempLSA_FoxHac09.csv";
	string H_prof_filename = "/home/rodney/git/corona3d_2020/src/inputs/Mars/H_density_profile_LSA_FoxHac09.csv";
	string Hplus_prof_filename = "/home/rodney/git/corona3d_2020/src/inputs/Mars/H+_density_profile_LSA_FoxHac09.csv";

	common::import_csv(temp_prof_filename, temp_profile[0], temp_profile[1], temp_profile[2], temp_profile[3]);
	common::import_csv(H_prof_filename, H_profile[0], H_profile[1]);
	common::import_csv(Hplus_prof_filename, Hplus_profile[0], Hplus_profile[1]);

	vector<double> logtemp;
	logtemp.resize(temp_profile[2].size());
	for (int i=0; i<logtemp.size(); i++)
	{
		logtemp[i] = log10(temp_profile[2][i]);
	}
	vector<double> logHplus;
	logHplus.resize(Hplus_profile[1].size());
	for (int i=0; i<logHplus.size(); i++)
	{
		logHplus[i] = log10(Hplus_profile[1][i]);
	}
	vector<double> logH;
	logH.resize(H_profile[1].size());
	for (int i=0; i<logH.size(); i++)
	{
		logH[i] = log10(H_profile[1][i]);
	}
	int num_alt_bins = 4000;

	vector<double> H_Hplus_rate;
	H_Hplus_rate.resize(num_alt_bins);
	alt_bins.resize(num_alt_bins);
	H_Hplus_CDF.resize(num_alt_bins);
	double rate_sum = 0.0;
	for (int i=0; i<num_alt_bins; i++)
	{
		alt_bins[i] = 2e7 + 10000.0*i;
		//double Ti = pow(10.0, common::interpolate(temp_profile[0], logtemp, alt_bins[i]));
		//double H_dens = pow(10.0, common::interpolate(H_profile[0], logH, alt_bins[i]));
		//double Hplus_dens = pow(10.0, common::interpolate(Hplus_profile[0], logHplus, alt_bins[i]));
		double Ti = common::interpolate(temp_profile[0], temp_profile[2], alt_bins[i]);
		double H_dens = common::interpolate(H_profile[0], H_profile[1], alt_bins[i]);
		double Hplus_dens = common::interpolate(Hplus_profile[0], Hplus_profile[1], alt_bins[i]);
		H_Hplus_rate[i] = H_Hplus_rate_coeff * sqrt(Ti) * H_dens * Hplus_dens;
		cout << alt_bins[i] << "," << H_Hplus_rate[i] << "\n";
		rate_sum = rate_sum + H_Hplus_rate[i];
	}
	for (int i=0; i<num_alt_bins; i++)
	{
		if (i == 0)
		{
			H_Hplus_CDF[i] = H_Hplus_rate[i] / rate_sum;
		}
		else
		{
			H_Hplus_CDF[i] = (H_Hplus_rate[i] / rate_sum) + H_Hplus_CDF[i-1];
		}
	}
}

Distribution_Hot_H::~Distribution_Hot_H() {

}

void Distribution_Hot_H::init(Particle* p)
{
	double my_mass = p->get_mass();

	// altitude distribution for hot H
	double r = get_new_radius_H_Hplus();
	double alt = r - my_planet.get_radius();
	double temp_ion = common::interpolate(temp_profile[0], temp_profile[2], alt);
	double temp_neut = common::interpolate(temp_profile[0], temp_profile[1], alt);

	double phi = constants::twopi*(common::get_rand());
	double u = 2.0*common::get_rand() - 1.0;
	double x = r*sqrt(1-(u*u))*cos(phi);
	double y = r*sqrt(1-(u*u))*sin(phi);
	double z = r*u;

	/*
	// Hemispherical Adjustment For Dayside Photochemical Process
	if (x < 0)
	{
		x = -x;
	}
	*/

	double vavg = sqrt(constants::k_b*temp_ion/m_Hion);  // average thermal H+ velocity
	double v_ion[] = {0.0, 0.0, 0.0};
	gen_mb(vavg, v_ion);


	/*
	// spherically isotropic velocity vector
	phi = constants::twopi*common.get_rand();
	u = 2.0*common.get_rand() - 1.0;
	double vx = v*sqrt(1-u*u)*cos(phi);
	double vy = v*sqrt(1-u*u)*sin(phi);
	double vz = v*u;

	// Add initial H+ and H translational momentum
	double vavg = sqrt(constants::k_b*temp_ion/(m_Hion));  // average thermal ion velocity
	double v_ion[] = {0.0, 0.0, 0.0};
	gen_mb(vavg, v_ion);                               // sample from Maxwell-Boltzmann distribution
	vavg = sqrt(constants::k_b*temp_neut/my_mass);    //average thermal H velocity
	double v_H[] = {0.0, 0.0, 0.0};
	gen_mb(vavg, v_H);                                 // sample from Maxwell-Boltzmann distribution

	// add it all up

	// case 1
	double vx = (m_Hion*v_ion[0] + my_mass*v_H[0]) / (m_Hion+my_mass);
	double vy = (m_Hion*v_ion[1] + my_mass*v_H[1]) / (m_Hion+my_mass);
	double vz = (m_Hion*v_ion[2] + my_mass*v_H[2]) / (m_Hion+my_mass);

	// case 2
	double vx = v_H[0] + (m_ion*v_ion[0] + my_mass*v_H[0]) / (m_ion+my_mass);
	double vy = v_H[1] + (m_ion*v_ion[1] + my_mass*v_H[1]) / (m_ion+my_mass);
	double vz = v_H[2] + (m_ion*v_ion[2] + my_mass*v_H[2]) / (m_ion+my_mass);
	*/

	p->init_particle(x, y, z, v_ion[0], v_ion[1], v_ion[2]);
}

// scans H_Hplus_CDF for new particle radius
double Distribution_Hot_H::get_new_radius_H_Hplus()
{
	double u = common::get_rand();
	int k = 0;
	while (H_Hplus_CDF[k] < u)
	{
		k++;
	}
	return alt_bins[k] + my_planet.get_radius();
}
