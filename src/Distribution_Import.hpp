/*
 * Distribution_Import.hpp
 *
 *  Created on: Jul 19, 2020
 *      Author: rodney
 */

#ifndef DISTRIBUTION_IMPORT_HPP_
#define DISTRIBUTION_IMPORT_HPP_

#include <iostream>
#include <sstream>
#include <fstream>
#include "Distribution.hpp"
using namespace std;

class Distribution_Import: public Distribution {
public:
	Distribution_Import();
	virtual ~Distribution_Import();
	void init(Particle* p);

private:
	int num_particles;
	int next_index;
	Matrix<double, Dynamic, 3> positions;
	Matrix<double, Dynamic, 3> velocities;
};

#endif /* DISTRIBUTION_IMPORT_HPP_ */
