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

int main(int argc, char* argv[])
{
	srand(time(NULL));

	int n = 10000;
	int timesteps = 100;
	double dt = 25.0;
	ofstream outFile;

	Particle testParticles[n];
	Planet venus(4.8675e24, 6051.8e3);
	initParticles(venus, testParticles, n);

	for (int i=0; i<timesteps; i++)
	{
		outFile.open("/home/rodney/Documents/coronaTest/data/" + to_string(i) + ".out");
		for (int j=0; j<n; j++)
		{
			outFile << setprecision(10) << testParticles[j].position[0] << '\t';
			outFile << setprecision(10) << testParticles[j].position[1] << '\t';
			outFile << setprecision(10) << testParticles[j].position[2] << '\n';
		}
		outFile.close();
		stepParticles(testParticles, dt, n, venus.k_g);
	}

	return 0;
}
