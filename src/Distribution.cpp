/*
 * Distribution.cpp
 *
 *  Created on: Jul 18, 2020
 *      Author: rodney
 */

#include "Distribution.hpp"

Distribution::Distribution(Planet my_p, double ref_h, double ref_T) {
	my_planet = my_p;
	ref_height = ref_h;
	ref_radius = my_planet.get_radius() + ref_height;
	ref_temp = ref_T;
}

Distribution::~Distribution() {

}
