/*
 * Distribution_Hot_O.cpp
 *
 *  Created on: Jul 22, 2020
 *      Author: rodney
 */

#include "Distribution_Hot_O.hpp"

Distribution_Hot_O::Distribution_Hot_O(Planet my_p, double ref_h, double ref_T)
	: Distribution(my_p, ref_h, ref_T) {

	DR160 = 600e12;
	H_DR = 19e5;
	T_ion = 400.0;
	T_e = 1600.0;
	m_ion = 31.9983*constants::amu;
}

Distribution_Hot_O::~Distribution_Hot_O() {

}

void Distribution_Hot_O::init(shared_ptr<Particle> p)
{
	// altitude distribution for O2+ dissociative recombination
	double r = my_planet.get_radius() + 160e5 - log(common::get_rand())*H_DR;

	double phi = constants::twopi*(common::get_rand());
	double u = 2.0*common::get_rand() - 1.0;
	double x = r*sqrt(1-(u*u))*cos(phi);
	double y = r*sqrt(1-(u*u))*sin(phi);
	double z = r*u;

	// Hemispherical Adjustment For Dayside Photochemical Process
	if (x < 0)
	{
		x = -x;
	}

	// Select Hot O Electronic Channel ###
	// Probabilities and energies
	// Particles are activated per channel to allow testing the contribution of each population in the simulation
	double Ei = 0.0;
	double randnum = common::get_rand();
	if (randnum < 0.22)
	{
		Ei = 6.99*constants::ergev;
	}
	else if (randnum < (0.22+0.42))
	{
		Ei = 5.02*constants::ergev;
	}
	else if (randnum < (0.22+0.42+0.31))
	{
		Ei = 3.06*constants::ergev;
	}
	else
	{
		Ei = 0.84*constants::ergev;
	}

	// Thermal Rotational Energy of O2+
	Ei = Ei + E_rot(3.35967e-16, T_ion);

	// Translational Energy per O Resulting from Dissociative Recombination of O2+
	double v = sqrt(Ei / p->get_mass());

	// spherically isotropic velocity vector
	phi = constants::twopi*common::get_rand();
	u = 2.0*common::get_rand() - 1.0;
	double vx = v*sqrt(1-u*u)*cos(phi);
	double vy = v*sqrt(1-u*u)*sin(phi);
	double vz = v*u;

	// Add initial ion and electron translational momentum
    double vavg = sqrt(constants::k_b*T_ion/(m_ion));  // average thermal ion velocity
    double v_ion[] = {0.0, 0.0, 0.0};
    gen_mb(vavg, v_ion);                               // sample from Maxwell-Boltzmann distribution
    vavg = sqrt(constants::k_b*T_e/constants::m_e);    //average thermal electron velocity
    double v_e[] = {0.0, 0.0, 0.0};
    gen_mb(vavg, v_e);                                 // sample from Maxwell-Boltzmann distribution

    // add it all up
    vx = vx + (m_ion*v_ion[0] + constants::m_e*v_e[0]) / (m_ion+constants::m_e);
    vy = vy + (m_ion*v_ion[1] + constants::m_e*v_e[1]) / (m_ion+constants::m_e);
    vz = vz + (m_ion*v_ion[2] + constants::m_e*v_e[2]) / (m_ion+constants::m_e);

    p->init_particle(x, y, z, vx, vy, vz);
}

// Sample Rotational Energy from a Thermal Population of Rigid Rotators
double Distribution_Hot_O::E_rot(double B, double T)
{
	int j = 0;
	int nj = 200;  // highest rotational level considered

	// partition function using two terms of Euler-MacLaurin summation formula
	double Q_rot = constants::k_b*T/B + 1.0/3.0;

	double P_tot = 0.0;
	double P_j = 0.0;
	double u = common::get_rand();

	for (int i=0; i<(nj+1); i++)
	{
		P_j = (2.0*j+1.0)*exp(-j*(j+1)*B/(constants::k_b*T))/Q_rot;
		P_tot = P_tot + P_j;
		j = i;
		if (u < P_tot)
		{
			break;
		}
	}
	return j*(j+1.0)*B;
}
