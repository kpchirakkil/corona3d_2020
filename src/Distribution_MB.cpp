/*
 * Distribution_MB.cpp
 *
 *  Created on: Jul 21, 2020
 *      Author: rodney
 */

#include "Distribution_MB.hpp"

Distribution_MB::Distribution_MB(Planet my_p, double ref_h, double ref_T)
	: Distribution(my_p, ref_h, ref_T) {

}

Distribution_MB::~Distribution_MB() {

}

void Distribution_MB::init(Particle* p)
{
	double phi = constants::twopi*(get_rand());
	double u = 2.0*get_rand() - 1;
	double v_avg = sqrt(constants::k_b*ref_temp/p->get_mass());

	double x = ref_radius*sqrt(1-(u*u))*cos(phi);
	double y = ref_radius*sqrt(1-(u*u))*sin(phi);
	double z = ref_radius*u;

	double randnum1 = get_rand();
	double randnum2 = get_rand();
	double randnum3 = get_rand();
	double randnum4 = get_rand();

	randnum1 = v_avg*sqrt(-2.0*log(1.0-randnum1));
	randnum2 = constants::twopi*randnum2;
	randnum3 = v_avg*sqrt(-2.0*log(1.0-randnum3));
	randnum4 = constants::twopi*randnum4;

	double vx = randnum1*cos(randnum2);
	double vy = randnum1*sin(randnum2);
	double vz = randnum3*cos(randnum4);

	p->init_particle(x, y, z, vx, vy, vz);
}

void Distribution_MB::init_vonly(Particle* p, double v_avg)
{
	double randnum1 = get_rand();
	double randnum2 = get_rand();
	double randnum3 = get_rand();
	double randnum4 = get_rand();

	randnum1 = v_avg*sqrt(-2.0*log(1.0-randnum1));
	randnum2 = constants::twopi*randnum2;
	randnum3 = v_avg*sqrt(-2.0*log(1.0-randnum3));
	randnum4 = constants::twopi*randnum4;

	double vx = randnum1*cos(randnum2);
	double vy = randnum1*sin(randnum2);
	double vz = randnum3*cos(randnum4);

	p->init_particle_vonly(vx, vy, vz);
}
