/*
 * Planet.hpp
 *
 *  Created on: May 29, 2020
 *      Author: rodney
 */

#ifndef PLANET_HPP_
#define PLANET_HPP_

class Planet {
public:
	Planet();
	void init(); // initialize with defaults (Venus mass/radius)
	void init(double m, double r);

	double mass;    // [kg] mass of planet
	double radius;  // [m] radius of planet
	double k_g;     // [m^3/s^2] planet's gravitational constant (-G*mass)
};

#endif /* PLANET_HPP_ */
