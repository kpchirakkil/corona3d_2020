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
#include <iostream>
#include <stdexcept>

Particle::Particle()
{
	active = true;
	traced = false;
	radius = 0.0;
	inverse_radius = 0.0;
	previous_radius = 0.0;
	position[0] = position[1] = position[2] = 0.0;
	velocity[0] = velocity[1] = velocity[2] = 0.0;
}

Particle::~Particle()
{

}

// deactivate this particle
void Particle::deactivate(string fate)
{
	active = false;

	// record fate of particle at bottom of collision log
	if (traced)
	{
		collision_log.push_back(fate);
	}
}

// All overloads use the same center-of-mass collision kernel.
void Particle::do_collision(shared_ptr<Particle> target, double theta, double time, double planet_r)
{
    do_collision(target, theta, time, planet_r, false, 0.0, -1, -1);
}

void Particle::do_collision(shared_ptr<Particle> target, double theta, double time, double planet_r, bool is_inelastic, double delta_E_eV)
{
    do_collision(target, theta, time, planet_r, is_inelastic, delta_E_eV, -1, -1);
}

void Particle::do_collision(shared_ptr<Particle> target, double theta, double time, double planet_r, bool is_inelastic, double delta_E_eV, int ji, int jf)
{
    if (!std::isfinite(theta) || theta<0 || theta>constants::pi || !std::isfinite(delta_E_eV))
        throw std::invalid_argument("Invalid collision angle or internal energy transfer");
    const double m=get_mass(), M=target->get_mass(), mu=m*M/(m+M);
    const Vector3d w(target->get_vx(),target->get_vy(),target->get_vz());
    const Vector3d g=velocity-w;
    const Vector3d cm=(m*velocity+M*w)/(m+M);
    const double gnorm=g.norm();
    const double energy=0.5*mu*g.squaredNorm()/constants::ergev;
    const double de=is_inelastic?delta_E_eV:0.0;
    // Closed channels must be removed by the sampler, never energy-clipped.
    if (de>energy) throw std::domain_error("Energetically closed inelastic transition");
    const double new_speed=std::sqrt(2.0*(energy-de)*constants::ergev/mu);
    Vector3d axis;
    if (gnorm>0) axis=g/gnorm;
    else {
        if (new_speed==0) return;
        // Exothermic zero-relative-speed limit has no preferred axis.
        const double z=2*common::get_rand()-1, phi=constants::twopi*common::get_rand();
        axis=Vector3d(std::sqrt(1-z*z)*std::cos(phi),std::sqrt(1-z*z)*std::sin(phi),z);
    }
    Vector3d seed=(std::abs(axis[2])<0.9)?Vector3d(0,0,1):Vector3d(1,0,0);
    Vector3d e1=(seed-axis*axis.dot(seed)).normalized();
    Vector3d e2(axis[1]*e1[2]-axis[2]*e1[1],axis[2]*e1[0]-axis[0]*e1[2],axis[0]*e1[1]-axis[1]*e1[0]);
    const double azimuth=constants::twopi*common::get_rand();
    const Vector3d direction=std::cos(theta)*axis+std::sin(theta)*(std::cos(azimuth)*e1+std::sin(azimuth)*e2);
    const double before=get_total_v()*1e-5;
    velocity=cm+(M/(m+M))*new_speed*direction;
    if (traced) {
        collision_log.push_back(to_string(time)+"\t\t"+to_string((radius-planet_r)*1e-5)+"\t"+target->get_name()+"\t"+
            to_string(theta*180/constants::pi)+"\t"+to_string(before)+"\t"+to_string(get_total_v()*1e-5)+"\t"+
            (is_inelastic?"INEL":"ELAS")+"\t"+to_string(de)+"\t"+to_string(ji)+"\t"+to_string(jf));
    }
}

void Particle::do_timestep(double dt, double k_g)
{
	previous_radius = radius;  // record current radius as new previous radius
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

// write collision log to given file
void Particle::dump_collision_log(string filename)
{
	ofstream outfile;
	outfile.open(filename);
	outfile << "#time(s)" << "\t\t" << "alt(km)" << "\t" << "targ" << "\t" << "angle(deg)" << "\t"
	        << "v_bef(km/s)" << "\t" << "v_aft(km/s)" << "\t" << "type" << "\t"
	        << "delta_E(eV)" << "\t" << "ji" << "\t" << "jf\n";
	int num_lines = collision_log.size();
	for (int i=0; i<num_lines; i++)
	{
		outfile << collision_log[i] << "\n";
	}
	outfile.close();
}

bool Particle::is_active() const
{
	return active;
}

bool Particle::is_traced() const
{
	return traced;
}

// return cosine of angle between particle trajectory and normal
double Particle::get_cos_theta() const
{
	double cos_theta = get_radial_v() / get_total_v();

	if (cos_theta > 1.0)
	{
		cos_theta = 1.0;
	}
	else if (cos_theta < -1.0)
	{
		cos_theta = -1.0;
	}

	return cos_theta;
}

double Particle::get_energy_in_eV() const
{
	return 0.5*get_mass()*pow(get_total_v(), 2.0)/constants::ergev;
}

double Particle::get_radial_energy_in_eV(double dt) const
{
	double radial_v = abs(radius - previous_radius) / dt;
	return 0.5*get_mass()*pow(radial_v, 2.0)/constants::ergev;
}

double Particle::get_radius() const
{
	return radius;
}

double Particle::get_radial_v() const
{
	double radial_v = (velocity[0]*position[0] + velocity[1]*position[1] + velocity[2]*position[2]) / radius;
	return radial_v;
}

double Particle::get_inverse_radius() const
{
	return inverse_radius;
}

double Particle::get_previous_radius() const
{
	return previous_radius;
}

double Particle::get_x() const
{
	return position[0];
}

double Particle::get_y() const
{
	return position[1];
}

double Particle::get_z() const
{
	return position[2];
}

double Particle::get_vx() const
{
	return velocity[0];
}

double Particle::get_vy() const
{
	return velocity[1];
}

double Particle::get_vz() const
{
	return velocity[2];
}

double Particle::get_total_v() const
{
	return sqrt(velocity[0]*velocity[0] +
			    velocity[1]*velocity[1] +
				velocity[2]*velocity[2]);
}

// initialize particle using given position and velocity
void Particle::init_particle(double x, double y, double z, double vx, double vy, double vz)
{
	radius = sqrt(x*x + y*y + z*z);
	inverse_radius = 1.0/radius;
	previous_radius = radius;
	position[0] = x;
	position[1] = y;
	position[2] = z;
	velocity[0] = vx;
	velocity[1] = vy;
	velocity[2] = vz;
}

// initialize particle with velocity only to be used for collisions
void Particle::init_particle_vonly(double vx, double vy, double vz)
{
	velocity[0] = vx;
	velocity[1] = vy;
	velocity[2] = vz;
}

// initialize a single particle at given radius using Maxwell-Boltzmann avg v
void Particle::init_particle_MB(double r, double v_avg)
{
	radius = r;
	double phi = constants::twopi*(common::get_rand());
	double u = 2.0*common::get_rand() - 1;
	inverse_radius = 1.0/r;
	previous_radius = radius;
	position[0] = r*sqrt(1-(u*u))*cos(phi);
	position[1] = r*sqrt(1-(u*u))*sin(phi);
	position[2] = r*u;

	double randnum1 = common::get_rand();
	double randnum2 = common::get_rand();
	double randnum3 = common::get_rand();
	double randnum4 = common::get_rand();

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
	double randnum1 = common::get_rand();
	double randnum2 = common::get_rand();
	double randnum3 = common::get_rand();
	double randnum4 = common::get_rand();

	randnum1 = v_avg*sqrt(-2.0*log(1.0-randnum1));
	randnum2 = constants::twopi*randnum2;
	randnum3 = v_avg*sqrt(-2.0*log(1.0-randnum3));
	randnum4 = constants::twopi*randnum4;

	velocity[0] = randnum1*cos(randnum2);
	velocity[1] = randnum1*sin(randnum2);
	velocity[2] = randnum3*cos(randnum4);
}

void Particle::set_traced()
{
	traced = true;
}
