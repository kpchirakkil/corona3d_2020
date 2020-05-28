/*
 * main.cpp
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */
#include <iostream>
#include "Particle.h"
#include "constants.h"

int main(int argc, char* argv[])
{
	Particle testParticle;
	testParticle.position[2] = jev;
	std::cout << testParticle.active << '\n';
	std::cout << testParticle.position[2] << '\n';

	return 0;
}
