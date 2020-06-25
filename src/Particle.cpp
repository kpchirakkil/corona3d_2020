/*
 * Particle.cpp
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */

#include "Particle.hpp"
#include "Particle_CO.hpp"
#include "Particle_CO2.hpp"
#include "Particle_H.hpp"
#include "Particle_N2.hpp"
#include "Particle_O.hpp"
#include "constants.hpp"

Particle::Particle()
{
	active = true;
	radius = 0.0;
	inverse_radius = 0.0;
	position[0] = position[1] = position[2] = 0.0;
	velocity[0] = velocity[1] = velocity[2] = 0.0;
}

Particle::~Particle()
{

}

// perform collision on a particle and update velocity vector
void Particle::do_collision(std::string targ_type, double targ_temp, double theta)
{
	double my_mass = get_mass();
	double targ_mass;
	Matrix<double, 3, 1> targ_v;
	Matrix<double, 3, 1> vcm;
	Matrix<double, 3, 3> Rrg;

	if (targ_type == "CO")
	{
		Particle_CO target;
		targ_mass = target.get_mass();
		target.init_particle_vonly_MB(sqrt(constants::k_b*targ_temp/targ_mass));
		targ_v = target.velocity;
	}
	else if (targ_type == "CO2")
	{
	   	Particle_CO2 target;
	   	targ_mass = target.get_mass();
	   	target.init_particle_vonly_MB(sqrt(constants::k_b*targ_temp/targ_mass));
	   	targ_v = target.velocity;
	}
	else if (targ_type == "H")
	{
		Particle_H target;
	    targ_mass = target.get_mass();
	    target.init_particle_vonly_MB(sqrt(constants::k_b*targ_temp/targ_mass));
	    targ_v = target.velocity;
	}
	else if (targ_type == "N2")
	{
	   	Particle_N2 target;
	   	targ_mass = target.get_mass();
	   	target.init_particle_vonly_MB(sqrt(constants::k_b*targ_temp/targ_mass));
	   	targ_v = target.velocity;
	}
	else if (targ_type == "O")
	{
	   	Particle_O target;
	   	targ_mass = target.get_mass();
	   	target.init_particle_vonly_MB(sqrt(constants::k_b*targ_temp/targ_mass));
	   	targ_v = target.velocity;
	}

	vcm = (my_mass*velocity.array() + targ_mass * targ_v.array()) / (my_mass + targ_mass);
	Matrix<double, 3, 1> v1v = velocity.array() - vcm.array();        // particle 1 c-o-m velocity
	Matrix<double, 3, 1> v2v = targ_v.array() - vcm.array();          // particle 2 c-o-m velocity
	double v1 = sqrt(v1v[0]*v1v[0] + v1v[1]*v1v[1] + v1v[2]*v1v[2]);  // particle 1 c-o-m scalar velocity
	double v2 = sqrt(v2v[0]*v2v[0] + v2v[1]*v2v[1] + v2v[2]*v2v[2]);  // particle 2 c-o-m scalar velocity

	// unit vector parallel to particle 1 velocity
	Matrix<double, 3, 1> r = velocity.array() / sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);

	double alpha = atan2(velocity[1], velocity[0]);
	double phi = atan2(velocity[2], sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1]));
	double gamma = constants::twopi*get_rand();

	Matrix<double, 3, 1> vp;
	vp[0] = v1*cos(alpha)*cos(phi-theta);
	vp[1] = v1*sin(alpha)*cos(phi-theta);
	vp[2] = v1*sin(phi-theta);

	double Cg = cos(gamma);
	double Sg = sin(gamma);
	double Vg = 1.0-Cg;


}

void Particle::do_timestep(double dt, double k_g)
{
	Array<double, 3, 1> a = {0.0, 0.0, 0.0}; // particle acceleration vector

	// calculate acceleration at current position
	double inv_r_cube = inverse_radius*inverse_radius*inverse_radius;
	a = k_g*position*inv_r_cube;

	// calculate next position and update particle
	position.array() = position.array() + (velocity.array()*dt) + (0.5*a*dt*dt);
	radius = sqrt(position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
	inverse_radius = 1.0 / radius;

	// calculate acceleration at next position
	inv_r_cube = inverse_radius*inverse_radius*inverse_radius;
	a = k_g*position.array()*inv_r_cube;

	// calculate next velocity using acceleration at next position and update particle
	velocity.array() = velocity.array() + 0.5*a*dt;
}

bool Particle::get_active()
{
	return active;
}

double Particle::get_x()
{
	return position[0];
}

double Particle::get_y()
{
	return position[1];
}

double Particle::get_z()
{
	return position[2];
}

double Particle::get_vx()
{
	return velocity[0];
}

double Particle::get_vy()
{
	return velocity[1];
}

double Particle::get_vz()
{
	return velocity[2];
}

// initialize a single particle at given radius using Maxwell-Boltzmann avg v
void Particle::init_particle_MB(double r, double v_avg)
{
	radius = r;
	double phi = constants::twopi*(get_rand());
	double u = 2.0*get_rand() - 1;
	inverse_radius = 1.0/r;
	position[0] = r*sqrt(1-(u*u))*cos(phi);
	position[1] = r*sqrt(1-(u*u))*sin(phi);
	position[2] = r*u;

	double randnum1 = get_rand();
	double randnum2 = get_rand();
	double randnum3 = get_rand();
	double randnum4 = get_rand();

	randnum1 = v_avg*sqrt(-2.0*log(1.0-randnum1));
	randnum2 = constants::twopi*randnum2;
	randnum3 = v_avg*sqrt(-2.0*log(1.0-randnum3));
	randnum4 = constants::twopi*randnum4;

	velocity[0] = randnum1*cos(randnum2);
	velocity[1] = randnum1*sin(randnum2);
	velocity[2] = randnum3*cos(randnum4);
}

void Particle::init_particle_vonly_MB(double v_avg)
{
	double randnum1 = get_rand();
	double randnum2 = get_rand();
	double randnum3 = get_rand();
	double randnum4 = get_rand();

	randnum1 = v_avg*sqrt(-2.0*log(1.0-randnum1));
	randnum2 = constants::twopi*randnum2;
	randnum3 = v_avg*sqrt(-2.0*log(1.0-randnum3));
	randnum4 = constants::twopi*randnum4;

	velocity[0] = randnum1*cos(randnum2);
	velocity[1] = randnum1*sin(randnum2);
	velocity[2] = randnum3*cos(randnum4);
}
