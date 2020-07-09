/*
 * Atmosphere.cpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "Atmosphere.hpp"
#include "constants.hpp"
using namespace std;

// construct atmosphere using specific parameters
Atmosphere::Atmosphere(int n, Planet p, Background_Species bg, double T, double ref_h)
{
	srand((unsigned)time(NULL));
	N = n;                      // number of test particles to track
	my_planet = p;
	my_parts.resize(N);
	T_bg = T;                   // [K] background temp where simulation starts
	ref_height = ref_h;         // [m] altitude for bottom of model
	bg_species = bg;

	double particle_r = my_planet.get_radius() + ref_height;
	double particle_vavg = sqrt(constants::k_b*T_bg/my_parts[0].get_mass());
	for (int i=0; i<N; i++)
	{
		my_parts[i].init_particle_MB(particle_r, particle_vavg);
	}
}

Atmosphere::~Atmosphere() {

}

// writes 3-column output file of all current particle positions
// file is saved to location specified by datapath
void Atmosphere::output_positions(string datapath)
{
	ofstream outfile;
	outfile.open(datapath);
	for (int i=0; i<N; i++)
	{
		outfile << setprecision(10) << my_parts[i].get_x() << '\t';
		outfile << setprecision(10) << my_parts[i].get_y() << '\t';
		outfile << setprecision(10) << my_parts[i].get_z() << '\n';
	}
	outfile.close();
}

// writes single-column output file of velocity bin counts in myParticles
void Atmosphere::output_velocity_distro(double bin_width, int num_bins, string datapath)
{
	ofstream outfile;
	outfile.open(datapath);
	double v = 0.0;             // velocity magnitude [m/s]
	int nvb = 0;                // bin number
	int vbins[num_bins] = {0};  // array of velocity bin counts

	for (int i=0; i<N; i++)
	{
		v = my_parts[i].get_total_v();
		nvb = (int)(v / bin_width);
		vbins[nvb]++;
	}
	for (int i=0; i<num_bins; i++)
	{
		outfile << vbins[i] << '\n';
	}
	outfile.close();
}

// iterate equation of motion and check for collisions for each active particle being tracked
void Atmosphere::run_simulation(double dt, int num_steps)
{
	double k = my_planet.get_k_g();

	for (int i=0; i<num_steps; i++)
	{
		output_positions("/home/rodney/Documents/coronaTest/data/positions" + to_string(i) + ".out");

		for (int j=0; j<N; j++)
		{
			if (my_parts[j].get_active() == false)
			{
				continue;
			}
			else
			{
				my_parts[j].do_timestep(dt, k);
				if (my_parts[j].get_active() == false)
				{
					continue;
				}
				else if (bg_species.check_collision(my_parts[j].get_radius(), my_parts[j].get_total_v(), dt))
				{
					my_parts[j].do_collision(bg_species.get_collision_target(), bg_species.get_collision_theta());
				}
			}
		}
	}

	std::cout << "Number of collisions: " << bg_species.get_num_collisions() << endl;
}
