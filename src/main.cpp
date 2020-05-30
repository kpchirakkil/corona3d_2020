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

	int n = 2000;
	double dt = 500.0;
	ofstream outFile1;
	ofstream outFile2;
	ofstream outFile3;
	outFile1.open("/home/rodney/Documents/corona3d.out1");
	outFile2.open("/home/rodney/Documents/corona3d.out2");
	outFile3.open("/home/rodney/Documents/corona3d.out3");

	Particle testParticles[n];
	Planet venus(4.8675e24, 6051.8e3);
	initParticles(venus, testParticles, n);

	for (int i=0; i<n; i++)
	{
		outFile1 << setprecision(10) << testParticles[i].position[0] << '\t';
		outFile1 << setprecision(10) << testParticles[i].position[1] << '\t';
		outFile1 << setprecision(10) << testParticles[i].position[2] << '\n';
	}
	outFile1.close();

	stepParticles(testParticles, dt, n, venus.k_g);

	for (int i=0; i<n; i++)
	{
		outFile2 << setprecision(10) << testParticles[i].position[0] << '\t';
		outFile2 << setprecision(10) << testParticles[i].position[1] << '\t';
		outFile2 << setprecision(10) << testParticles[i].position[2] << '\n';
	}
	outFile2.close();

	stepParticles(testParticles, dt, n, venus.k_g);

	for (int i=0; i<n; i++)
	{
		outFile3 << setprecision(10) << testParticles[i].position[0] << '\t';
		outFile3 << setprecision(10) << testParticles[i].position[1] << '\t';
		outFile3 << setprecision(10) << testParticles[i].position[2] << '\n';
	}
	outFile3.close();

	return 0;
}
