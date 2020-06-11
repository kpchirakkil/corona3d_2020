/*
 * Atmosphere.hpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#ifndef ATMOSPHERE_HPP_
#define ATMOSPHERE_HPP_

#include "Planet.hpp"
#include "Particle.hpp"

class Atmosphere {
public:
	Atmosphere(); // initialize with defaults (10K MB-distributed H atoms 160km above Venus)
	Atmosphere(int n, double planet_mass, double planet_radius, double particle_m, double T, double model_b);
	void doCollision(int paticle_idx, double m2);
	void doTimestep(double dt);
	void output_positions(std::string datapath);
	void output_velocity_distro(double bin_width, int num_bins, std::string datapath);
	void runSimulation(double dt, int numSteps);

	int N;                  // number of particles to track
	Planet *myPlanet;       // contains planet mass and radius
	Particle *myParticles;  // array of particles to be tracked
	double particleMass;    // [kg] mass of particle being tracked
	double T_bg;            // [K] background temp where simulation starts
	double model_bottom;    // [m] altitude above planet surface of model bottom
};

#endif /* ATMOSPHERE_HPP_ */
