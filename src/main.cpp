/*
 * main.cpp
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "Planet.h"
#include "Particle.h"
#include "constants.h"
#include "phys.h"
using namespace std;

// writes 3-column output file of all current particle positions in parts[]
// file is saved to location specified by datapath
void output_positions(Particle parts[], int n, string datapath)
{
	ofstream outFile;
	outFile.open(datapath);
	for (int i=0; i<n; i++)
	{
		outFile << setprecision(10) << parts[i].position[0] << '\t';
		outFile << setprecision(10) << parts[i].position[1] << '\t';
		outFile << setprecision(10) << parts[i].position[2] << '\n';
	}
	outFile.close();
}

void output_velocity_distro(Particle parts[], int n, string datapath)
{
	ofstream outFile;
	outFile.open(datapath);
	double v = 0.0;
	double dv = 100.0;
	int nvb = 0;
	int nbins = 100;
	int vbins[nbins] = {0};

	for (int i=0; i<n; i++)
	{
		v = sqrt(parts[i].velocity[0]*parts[i].velocity[0] +
				 parts[i].velocity[1]*parts[i].velocity[1] +
				 parts[i].velocity[2]*parts[i].velocity[2]);
		nvb = (int)(v / dv);
		vbins[nvb]++;
	}
	for (int i=0; i<nbins; i++)
	{
		outFile << vbins[i] << '\n';
	}
	outFile.close();
}

int main(int argc, char* argv[])
{
	srand(time(NULL));

	int n = 10000;
	int timesteps = 1000;
	double dt = 0.5;

	Particle testParticles[n];
	Planet venus(4.8675e24, 6051.8e3, 350.0, 1.00794*amu);
	initParticles(venus, testParticles, n);

	output_velocity_distro(testParticles, n, "/home/rodney/Documents/coronaTest/vdist.out");

	for (int i=0; i<timesteps; i++)
	{

		output_positions(testParticles, n, "/home/rodney/Documents/coronaTest/data/positions" + to_string(i) + ".out");
		stepParticles(testParticles, dt, n, venus.k_g);
	}

	output_velocity_distro(testParticles, n, "/home/rodney/Documents/coronaTest/vdist2.out");

	return 0;
}
