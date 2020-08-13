/*
 * Distribution_Hot_H.cpp
 *
 *  Created on: Jul 31, 2020
 *      Author: rodney
 */

#include "Distribution_Hot_H.hpp"

Distribution_Hot_H::Distribution_Hot_H(Planet my_p, double ref_h, double ref_T)
	: Distribution(my_p, ref_h, ref_T) {

	T_ion = 400.0;
	T_e = 1250.0;
	m_ion = 1.00728*constants::amu;
}

Distribution_Hot_H::~Distribution_Hot_H() {

}

void Distribution_Hot_H::init(Particle* p)
{
	double my_mass = p->get_mass();

	// altitude distribution for hot H
	double r = my_planet.get_radius() + ref_height;

	double phi = constants::twopi*(get_rand());
	double u = 2.0*get_rand() - 1.0;
	double x = r*sqrt(1-(u*u))*cos(phi);
	double y = r*sqrt(1-(u*u))*sin(phi);
	double z = r*u;

	/*
	// Hemispherical Adjustment For Dayside Photochemical Process
	if (x < 0)
	{
		x = -x;
	}

	// spherically isotropic velocity vector
	phi = constants::twopi*get_rand();
	u = 2.0*get_rand() - 1.0;
	double vx = v*sqrt(1-u*u)*cos(phi);
	double vy = v*sqrt(1-u*u)*sin(phi);
	double vz = v*u;
	*/

	// Add initial H+ and H translational momentum
	double vavg = sqrt(constants::k_b*T_ion/(m_ion));  // average thermal ion velocity
	double v_ion[] = {0.0, 0.0, 0.0};
	gen_mb(vavg, v_ion);                               // sample from Maxwell-Boltzmann distribution
	vavg = sqrt(constants::k_b*ref_temp/my_mass);    //average thermal H velocity
	double v_H[] = {0.0, 0.0, 0.0};
	gen_mb(vavg, v_H);                                 // sample from Maxwell-Boltzmann distribution

	// add it all up

	// case 1
	double vx = (m_ion*v_ion[0] + my_mass*v_H[0]) / (m_ion+my_mass);
	double vy = (m_ion*v_ion[1] + my_mass*v_H[1]) / (m_ion+my_mass);
	double vz = (m_ion*v_ion[2] + my_mass*v_H[2]) / (m_ion+my_mass);

	/* case 2
	double vx = v_H[0] + (m_ion*v_ion[0] + my_mass*v_H[0]) / (m_ion+my_mass);
	double vy = v_H[1] + (m_ion*v_ion[1] + my_mass*v_H[1]) / (m_ion+my_mass);
	double vz = v_H[2] + (m_ion*v_ion[2] + my_mass*v_H[2]) / (m_ion+my_mass);
	*/

	p->init_particle(x, y, z, vx, vy, vz);
}

