/*
 * Particle.cpp
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */

#include "Particle.hpp"
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
