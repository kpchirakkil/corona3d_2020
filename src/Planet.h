/*
 * Planet.h
 *
 *  Created on: May 29, 2020
 *      Author: rodney
 */

#ifndef PLANET_H_
#define PLANET_H_

class Planet {
public:
	Planet(double m, double r, double T, double part_m);
	double mass;        // [kg] mass of planet
	double radius;      // [m] radius of planet
	double k_g;         // [m^3/s^2] planet's gravitational constant (-G*mass)
	double particle_m;  // [kg] mass of particle to be tracked
	double T_bg;        // [K] background temp for MB generation
};

#endif /* PLANET_H_ */
