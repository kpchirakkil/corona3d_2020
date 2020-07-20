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
#include <iostream>

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

// deactivate this particle
void Particle::deactivate()
{
	active = false;
}

// perform collision on a particle and update velocity vector
void Particle::do_collision(Particle* target, double theta)
{
	double my_mass = get_mass();
	double targ_mass = target->get_mass();
	Matrix<double, 3, 1> targ_v = {target->get_vx(), target->get_vy(), target->get_vz()};
	Matrix<double, 3, 1> vcm;
	Matrix<double, 3, 3> Rrg;

	vcm = (my_mass*velocity.array() + targ_mass * targ_v.array()) / (my_mass + targ_mass);
	Matrix<double, 3, 1> v1v = velocity.array() - vcm.array();        // particle 1 c-o-m velocity
	// Matrix<double, 3, 1> v2v = targ_v.array() - vcm.array();          // particle 2 c-o-m velocity
	double v1 = sqrt(v1v[0]*v1v[0] + v1v[1]*v1v[1] + v1v[2]*v1v[2]);  // particle 1 c-o-m scalar velocity
	// double v2 = sqrt(v2v[0]*v2v[0] + v2v[1]*v2v[1] + v2v[2]*v2v[2]);  // particle 2 c-o-m scalar velocity

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

	Rrg(0, 0) = r[0]*r[0]*Vg+Cg;
	Rrg(0, 1) = r[0]*r[1]*Vg+r[2]*Sg;
	Rrg(0, 2) = r[0]*r[2]*Vg-r[1]*Sg;

	Rrg(1, 0) = r[0]*r[1]*Vg-r[2]*Sg;
	Rrg(1, 1) = r[1]*r[1]*Vg+Cg;
	Rrg(1, 2) = r[1]*r[2]*Vg+r[0]*Sg;

	Rrg(2, 0) = r[0]*r[2]*Vg+r[1]*Sg;
	Rrg(2, 1) = r[1]*r[2]*Vg-r[0]*Sg;
	Rrg(2, 2) = r[2]*r[2]*Vg+Cg;

	Matrix<double, 3, 1> vrel1 = Rrg * vp;

	// update post-collision velocity
	velocity = vcm.array() + vrel1.array();

	// in case you need the updated collision partner velocity for something
	// targ_v = vcm.array() - ((my_mass / targ_mass) * vrel1.array());
}

void Particle::do_timestep(double dt, double k_g)
{
	Array<double, 3, 1> a = {0.0, 0.0, 0.0}; // particle acceleration vector

	// calculate acceleration at current position
	double inv_r_cube = inverse_radius*inverse_radius*inverse_radius;
	a = k_g*position.array()*inv_r_cube;

	// calculate next position and update particle
	position.array() = position.array() + (velocity.array()*dt) + (0.5*a*dt*dt);
	radius = sqrt(position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
	inverse_radius = 1.0 / radius;

	// calculate acceleration at next position
	inv_r_cube = inverse_radius*inverse_radius*inverse_radius;
	a = a + k_g*position.array()*inv_r_cube;

	// calculate next velocity using acceleration at next position and update particle
	velocity.array() = velocity.array() + 0.5*a*dt;
}

bool Particle::get_active()
{
	return active;
}

double Particle::get_radius()
{
	return radius;
}

double Particle::get_inverse_radius()
{
	return inverse_radius;
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

double Particle::get_total_v()
{
	return sqrt(velocity[0]*velocity[0] +
			    velocity[1]*velocity[1] +
				velocity[2]*velocity[2]);
}

// initialize particle using custom position and velocity
void Particle::init_particle(double x, double y, double z, double vx, double vy, double vz)
{
	radius = sqrt(x*x + y*y + z*z);
	inverse_radius = 1.0/radius;
	position[0] = x;
	position[1] = y;
	position[2] = z;
	velocity[0] = vx;
	velocity[1] = vy;
	velocity[2] = vz;
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
