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
	inline const double pi    = M_PIl;            // pi [unitless]
	inline const double twopi = 2*pi;             // 2*pi [unitless]
	inline const double k_b   = 1.3806504e-23;    // Boltzmann's Constant [J/K]
	inline const double c     = 299792458.0;      // Speed of Light in Vacuum [m/s]
	inline const double G     = 6.67428e-11;      // Gravitational Constant [m^3/kg/s^2]
	inline const double amu   = 1.660538782e-27;  // Atomic Mass Unit [kg]
	inline const double m_e   = 9.10938215e-31;   // Electron Mass [kg]
	inline const double q_e   = 1.602176487e-19;  // Elementary Charge [C]
	inline const double jev   = q_e;              // Joules/Electron Volt [unitless]
}

#endif /* CONSTANTS_HPP_ */
