/*
 * Atmosphere.cpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "Atmosphere.hpp"
#include "constants.hpp"
using namespace std;

// default constructor (makes 10K MB-distributed H atoms 160km above Venus)
Atmosphere::Atmosphere()
{
	srand(time(NULL));
	N = 10000;
	myPlanet = new Planet;
	myPlanet->init();               // Venus mass and radius by default
	myParticles = new Particle[N];  // use default of 10000 particles
	particleMass = 1.00794*amu;     // [kg] mass of hydrogen atom
	T_bg = 350.0;                   // [K] background temp where simulation starts
	model_bottom = 160e3;           // [m] default altitude for bottom of model

	double particle_r = myPlanet->radius + model_bottom;
	double particle_vavg = sqrt(k_b*T_bg/particleMass);
	for (int i=0; i<N; i++)
	{
		myParticles[i].initParticle_MB(particle_r, particle_vavg);
	}
}

// construct atmosphere using specific parameters
Atmosphere::Atmosphere(int n, double planet_mass, double planet_radius, double particle_m, double T, double model_b)
{
	srand(time(NULL));
	N = n;
	myPlanet = new Planet;
	myPlanet->init(planet_mass, planet_radius);
	myParticles = new Particle[N];
	particleMass = particle_m;
	T_bg = T;
	model_bottom = model_b;

	double particle_r = myPlanet->radius + model_bottom;
	double particle_vavg = sqrt(k_b*T_bg/particleMass);
	for (int i=0; i<N; i++)
	{
		myParticles[i].initParticle_MB(particle_r, particle_vavg);
	}
}

// perform collision on a particle and update velocity vector
void Atmosphere::doCollision(int paticle_idx, double m2)
{

}

// iterate equation of motion for all particles using velocity Verlet algorithm
void Atmosphere::doTimestep(double dt)
{
	double a[3] = {0.0, 0.0, 0.0}; // particle acceleration vector
	for (int i=0; i<N; i++)
	{
		if (myParticles[i].active == false)
		{
			continue;
		}
		else
		{
			// calculate acceleration at current position
			double invRcube = pow(myParticles[i].inverseRadius, 3);
			a[0] = myPlanet->k_g*myParticles[i].position[0]*invRcube;
			a[1] = myPlanet->k_g*myParticles[i].position[1]*invRcube;
			a[2] = myPlanet->k_g*myParticles[i].position[2]*invRcube;

			// calculate next position and update particle
			myParticles[i].position[0] = myParticles[i].position[0] + myParticles[i].velocity[0]*dt + 0.5*a[0]*dt*dt;
			myParticles[i].position[1] = myParticles[i].position[1] + myParticles[i].velocity[1]*dt + 0.5*a[1]*dt*dt;
			myParticles[i].position[2] = myParticles[i].position[2] + myParticles[i].velocity[2]*dt + 0.5*a[2]*dt*dt;
			myParticles[i].radius = sqrt(pow(myParticles[i].position[0], 2) + pow(myParticles[i].position[1], 2) + pow(myParticles[i].position[2], 2));
			myParticles[i].inverseRadius = 1 / myParticles[i].radius;

			// calculate acceleration at next position
			invRcube = pow(myParticles[i].inverseRadius, 3);
			a[0] = myPlanet->k_g*myParticles[i].position[0]*invRcube;
			a[1] = myPlanet->k_g*myParticles[i].position[1]*invRcube;
			a[2] = myPlanet->k_g*myParticles[i].position[2]*invRcube;

			// calculate next velocity using acceleration at next position and update particle
			myParticles[i].velocity[0] = myParticles[i].velocity[0] + 0.5*a[0]*dt;
			myParticles[i].velocity[1] = myParticles[i].velocity[1] + 0.5*a[1]*dt;
			myParticles[i].velocity[2] = myParticles[i].velocity[2] + 0.5*a[2]*dt;
		}
	}
}

// writes 3-column output file of all current particle positions
// file is saved to location specified by datapath
void Atmosphere::output_positions(string datapath)
{
	ofstream outFile;
	outFile.open(datapath);
	for (int i=0; i<N; i++)
	{
		outFile << setprecision(10) << myParticles[i].position[0] << '\t';
		outFile << setprecision(10) << myParticles[i].position[1] << '\t';
		outFile << setprecision(10) << myParticles[i].position[2] << '\n';
	}
	outFile.close();
}

// writes single-column output file of velocity bin counts in myParticles
void Atmosphere::output_velocity_distro(double bin_width, int num_bins, string datapath)
{
	ofstream outFile;
	outFile.open(datapath);
	double v = 0.0;             // velocity magnitude [m/s]
	int nvb = 0;                // bin number
	int vbins[num_bins] = {0};  // array of velocity bin counts

	for (int i=0; i<N; i++)
	{
		v = sqrt(myParticles[i].velocity[0]*myParticles[i].velocity[0] +
				 myParticles[i].velocity[1]*myParticles[i].velocity[1] +
				 myParticles[i].velocity[2]*myParticles[i].velocity[2]);
		nvb = (int)(v / bin_width);
		vbins[nvb]++;
	}
	for (int i=0; i<num_bins; i++)
	{
		outFile << vbins[i] << '\n';
	}
	outFile.close();
}

void runSimulation(double dt, int numSteps)
{

}
