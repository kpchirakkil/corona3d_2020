/*
 * constants.hpp
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */

#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_

#include <cmath>

namespace constants
{
	const double pi    = M_PIl;            // pi [unitless]
	const double twopi = 2*pi;             // 2*pi [unitless]
	const double k_b   = 1.3806504e-16;    // Boltzmann's Constant [erg/K]
	const double c     = 29979245800.0;    // Speed of Light in Vacuum [cm/s]
	const double G     = 6.67430e-8;       // Gravitational Constant [cm^3/g/s^2]
	const double amu   = 1.660538782e-24;  // Atomic Mass Unit [g]
	const double m_e   = 9.10938215e-28;   // Electron Mass [g]
	const double q_e   = 1.602176487e-19;  // Elementary Charge [C]
	const double jev   = q_e;              // Joules/Electron Volt [unitless]
	const double ergev = jev*1.0e7;        // ergs/Electron Volt [unitless]
}

#endif /* CONSTANTS_HPP_ */
