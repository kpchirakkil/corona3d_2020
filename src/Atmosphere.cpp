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
	srand((unsigned)time(NULL));
	N = 10000;
	Planet myPlanet;
	myPlanet.init();
	myH.resize(N);
	T_bg = 350.0;                           // [K] background temp where simulation starts
	model_bottom = 160e3;                   // [m] default altitude for bottom of model

	double particle_r = myPlanet.get_radius() + model_bottom;
	double particle_vavg = sqrt(constants::k_b*T_bg/myH[0].mass);
	for (int i=0; i<N; i++)
	{
		myH[i].init_particle_MB(particle_r, particle_vavg);
	}
}

// construct atmosphere using specific parameters
Atmosphere::Atmosphere(int n, Planet p, double T, double model_b)
{
	srand((unsigned)time(NULL));
	N = n;
	myPlanet = p;
	myH.resize(N);
	T_bg = T;
	model_bottom = model_b;

	double particle_r = myPlanet.get_radius() + model_bottom;
	double particle_vavg = sqrt(constants::k_b*T_bg/myH[0].mass);
	for (int i=0; i<N; i++)
	{
		myH[i].init_particle_MB(particle_r, particle_vavg);
	}
}

Atmosphere::~Atmosphere() {

}

// iterate equation of motion for each active particle being tracked
void Atmosphere::do_timestep(double dt)
{
	double k = myPlanet.get_k_g();
	for (int i=0; i<N; i++)
	{
		if (myH[i].get_active() == false)
		{
			continue;
		}
		else
		{
			myH[i].do_timestep(dt, k);
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
		outFile << setprecision(10) << myH[i].get_x() << '\t';
		outFile << setprecision(10) << myH[i].get_y() << '\t';
		outFile << setprecision(10) << myH[i].get_z() << '\n';
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
		v = sqrt(myH[i].get_vx()*myH[i].get_vx() +
				 myH[i].get_vy()*myH[i].get_vy() +
				 myH[i].get_vz()*myH[i].get_vz());
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
