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
	ofstream outFile;
	outFile.open("/home/rodney/Documents/corona3d.out");
	Particle testParticles[n];
	Planet venus(4.8675e24, 6051.8e3);


	initParticles(venus, testParticles, n);
	for (int i=0; i<n; i++)
	{
		outFile << setprecision(10) << testParticles[i].position[0] << '\t';
		outFile << setprecision(10) << testParticles[i].position[1] << '\t';
		outFile << setprecision(10) << testParticles[i].position[2] << '\n';
	}
	outFile.close();

	return 0;
}
