/*
 * Common_Functions.hpp
 *
 *  Created on: Aug 21, 2020
 *      Author: rodney
 */

#ifndef COMMON_FUNCTIONS_HPP_
#define COMMON_FUNCTIONS_HPP_

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
using namespace std;

class Common_Functions {
public:
	Common_Functions();
	virtual ~Common_Functions();

	void import_csv(string filename, vector<double> &col1, vector<double> &col2);
	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3);
	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4);
	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4, vector<double> &col5);
	void import_csv(string filename, vector<double> &col1, vector<double> &col2, vector<double> &col3, vector<double> &col4, vector<double> &col5, vector<double> &col6);

	// returns interpolated value at x from parallel arrays (x_data, y_data)
	double interpolate(vector<double> &x_data, vector<double> &y_data, double x);

	// returns uniformly distributed random number between 0 and 1
	double get_rand();
};

#endif /* COMMON_FUNCTIONS_HPP_ */
