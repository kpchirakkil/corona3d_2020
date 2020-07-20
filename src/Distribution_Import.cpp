/*
 * Distribution_Import.cpp
 *
 *  Created on: Jul 19, 2020
 *      Author: rodney
 */

#include "Distribution_Import.hpp"

Distribution_Import::Distribution_Import() {
	num_particles = 0;
	next_index = 0;

	ifstream pos_file, vel_file;
	pos_file.open("/home/rodney/Documents/research/software/Deighan/corona3d_2020/positions.out");
	vel_file.open("/home/rodney/Documents/research/software/Deighan/corona3d_2020/velocities.out");
	string pos_line, vel_line;

	while (getline(pos_file, pos_line) && getline(vel_file, vel_line))
	{
		num_particles++;
		stringstream str_pos(pos_line);
		stringstream str_vel(vel_line);

		positions.conservativeResize(positions.rows()+1, NoChange);
		velocities.conservativeResize(velocities.rows()+1, NoChange);

		str_pos >> positions(positions.rows()-1, 0) >> positions(positions.rows()-1, 1) >> positions(positions.rows()-1, 2);
		str_vel >> velocities(velocities.rows()-1, 0) >> velocities(velocities.rows()-1, 1) >> velocities(velocities.rows()-1, 2);
	}
	pos_file.close();
	vel_file.close();
}

Distribution_Import::~Distribution_Import() {

}

void Distribution_Import::init(Particle* p)
{
	if (next_index < num_particles)
	{
		p->init_particle(positions(next_index, 0), positions(next_index, 1), positions(next_index, 2), velocities(next_index, 0), velocities(next_index, 1), velocities(next_index, 2));
		next_index++;
	}
	else
	{
		cout << "End of available particles in imported distribution met!" << endl;
	}
}
