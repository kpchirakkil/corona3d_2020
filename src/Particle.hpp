/*
 * Particle.hpp
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */

#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_

#include <cstdlib>
#include <eigen3/Eigen/Core>
using namespace Eigen;

class Particle {
public:
	Particle();
	virtual ~Particle();
	virtual double get_mass() = 0;  //must be implemented in derived classes

	void deactivate();
	void do_collision(Particle* target, double theta);
	void do_timestep(double dt, double k_g);
	bool get_active();
	double get_radius();
	double get_inverse_radius();
	double get_x();
	double get_y();
	double get_z();
	double get_vx();
	double get_vy();
	double get_vz();
	double get_total_v();
	void init_particle(double x, double y, double z, double vx, double vy, double vz);
	void init_particle_vonly(double vx, double vy, double vz);
	void init_particle_MB(double r, double v_avg); // init particle from MB distribution
	void init_particle_vonly_MB(double v_avg);     // init with velocity only for collision partners

protected:
	bool active;                    // flag for whether particle is active, i.e. should still be considered in the simulation
	double radius;                  // radius from center of planet	[m]
	double inverse_radius;          // inverse radius (for computational efficiency) [m^-1]
	Matrix<double, 3, 1> position;  // position vector [m,m,m]
	Matrix<double, 3, 1> velocity;  // velocity vector [m/s,m/s,m/s]

	// returns uniformly distributed random number between 0 and 1
	double get_rand() {return ((double)rand() / (double)RAND_MAX);}
};

#endif /* PARTICLE_HPP_ */
