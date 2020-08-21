/*
 * Common_Functions.hpp
 *
 *  Created on: Aug 21, 2020
 *      Author: rodney
 */

#ifndef COMMON_FUNCTIONS_HPP_
#define COMMON_FUNCTIONS_HPP_

#include <vector>
#include <cstdlib>
using namespace std;

class Common_Functions {
public:
	Common_Functions();
	virtual ~Common_Functions();

	// returns interpolated value at x from parallel arrays (x_data, y_data)
	double interpolate(vector<double> &x_data, vector<double> &y_data, double x);

	// returns uniformly distributed random number between 0 and 1
	double get_rand();
};

#endif /* COMMON_FUNCTIONS_HPP_ */
