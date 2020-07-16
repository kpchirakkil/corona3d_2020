/*
 * Atmosphere.cpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#include "Atmosphere.hpp"
using namespace std;

// construct atmosphere using specific parameters
Atmosphere::Atmosphere(int n, Planet p, Background_Species bg, double T, double ref_h)
{
	srand((unsigned)time(NULL));
	N = n;                      // number of test particles to track
	active_parts = N;
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

	// everything below should eventually be incorporated into a Distribution class
	/*
	ifstream pos_file, vel_file;
	pos_file.open("/home/rodney/Documents/research/software/Deighan/corona3d_2020/positions2.out");
	vel_file.open("/home/rodney/Documents/research/software/Deighan/corona3d_2020/velocities2.out");
	string pos_line, vel_line;
	double x, y, z, vx, vy, vz;

	for (int i=0; i<N; i++)
	{
		getline(pos_file, pos_line);
		getline(vel_file, vel_line);
		stringstream str_pos(pos_line);
		stringstream str_vel(vel_line);
		str_pos >> x >> y >> z;
		str_vel >> vx >> vy >> vz;
		my_parts[i].init_particle_custom(x, y, z, vx, vy, vz);
	}
	pos_file.close();
	vel_file.close();
	*/
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
		if (my_parts[i].get_active())
		{
			v = my_parts[i].get_total_v();
			nvb = (int)(v / bin_width);
			vbins[nvb]++;
		}
	}
	for (int i=0; i<num_bins; i++)
	{
		outfile << vbins[i] << '\n';
	}
	outfile.close();
}

// iterate equation of motion and check for collisions for each active particle being tracked
// a lot of stuff in here needs to be changed to be dynamically determined at runtime
void Atmosphere::run_simulation(double dt, int num_steps)
{
	double k = my_planet.get_k_g();
	double v_Obg = sqrt(8.0*constants::k_b*T_bg / (constants::pi*my_parts[0].get_mass()));
	//double v_Obg = 0;
	//double thresh_v = sqrt(2.0*constants::G*my_planet.get_mass()*(my_parts[0].get_inverse_radius()-1.0/(my_planet.get_radius()+900e3)));

	for (int i=0; i<num_steps; i++)
	{
		//output_positions("/home/rodney/Documents/coronaTest/data/positions" + to_string(i) + ".out");

		for (int j=0; j<N; j++)
		{
			if (my_parts[j].get_active())
			{
				my_parts[j].do_timestep(dt, k);
				if (bg_species.check_collision(my_parts[j].get_radius(), my_parts[j].get_total_v(), dt))
				{
					my_parts[j].do_collision(bg_species.get_collision_target(), bg_species.get_collision_theta());
				}
				if (my_parts[j].get_radius() < (my_planet.get_radius() + 900e3) && (my_parts[j].get_total_v() + v_Obg) < sqrt(2.0*constants::G*my_planet.get_mass()*(my_parts[j].get_inverse_radius()-1.0/(my_planet.get_radius()+900e3))))
				{
					my_parts[j].deactivate();
					active_parts--;
				}
			}
		}
	}

	std::cout << "Number of collisions: " << bg_species.get_num_collisions() << endl;
	std::cout << "Active particles remaining: " << active_parts << endl;
}
