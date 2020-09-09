/*
 * Atmosphere.cpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#include "Atmosphere.hpp"

// construct atmosphere using given parameters
Atmosphere::Atmosphere(int n, Planet p, vector<Particle*> parts, Distribution* dist, Background_Species bg, double T, double ref_h, string temp_profile)
{
	srand((unsigned)time(NULL));  // seed random number generator
	num_parts = n;                // number of test particles to track
	active_parts = num_parts;
	my_planet = p;
	my_dist = dist;
	my_parts.resize(num_parts);
	T_bg = T;                     // [K] background temp where simulation starts
	ref_height = ref_h;           // [cm] altitude for bottom of model
	bg_species = bg;

	for (int i=0; i<num_parts; i++)
	{
		my_parts[i] = parts[i];
		my_dist->init(my_parts[i]);
	}

	//the way I initialized particles before the distribution class was made
	//only works for MB dist, but leaving in just in case it comes in handy
	/*
	double particle_r = my_planet.get_radius() + ref_height;
	double particle_vavg = sqrt(constants::k_b*T_bg/my_parts[0].get_mass());
	for (int i=0; i<N; i++)
	{
		my_parts[i].init_particle_MB(particle_r, particle_vavg);
	}
	*/

	// read in temperature profile
	common.import_csv(temp_profile, alt_bins, Tn, Ti, Te);
}

Atmosphere::~Atmosphere() {

}

// writes 3-column output file of all current particle positions
// file is saved to location specified by datapath
void Atmosphere::output_positions(string datapath)
{
	ofstream outfile;
	outfile.open(datapath);
	for (int i=0; i<num_parts; i++)
	{
		outfile << setprecision(10) << my_parts[i]->get_x() << '\t';
		outfile << setprecision(10) << my_parts[i]->get_y() << '\t';
		outfile << setprecision(10) << my_parts[i]->get_z() << '\n';
	}
	outfile.close();
}

// writes single-column output file of altitude bin counts using active particles
// bin_width is in cm
void Atmosphere::output_altitude_distro(double bin_width, int num_bins, std::string datapath)
{
	ofstream outfile;
	outfile.open(datapath);
	double alt = 0.0;           // altitude of particle [cm]
	int nb = 0;                 // bin number
	int abins[num_bins] = {0};  // array of altitude bin counts

	for (int i=0; i<num_parts; i++)
	{
		if (my_parts[i]->get_active())
		{
			alt = my_parts[i]->get_radius() - my_planet.get_radius();
			nb = (int)(alt / bin_width);
			abins[nb]++;
		}
	}
	for (int i=0; i<num_bins; i++)
	{
		outfile << abins[i] << "\n";
	}
	outfile.close();
}

// writes single-column output file of velocity bin counts using active particles
// bin_width is in cm/s
void Atmosphere::output_velocity_distro(double bin_width, int num_bins, string datapath)
{
	ofstream outfile;
	outfile.open(datapath);
	double v = 0.0;             // velocity magnitude [cm/s]
	int nb = 0;                 // bin number
	int vbins[num_bins] = {0};  // array of velocity bin counts

	for (int i=0; i<num_parts; i++)
	{
		if (my_parts[i]->get_active())
		{
			v = my_parts[i]->get_total_v();
			nb = (int)(v / bin_width);
			vbins[nb]++;
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
	double v_Obg = sqrt(8.0*constants::k_b*T_bg / (constants::pi*15.9994*constants::amu));
	//double v_Obg = 0;
	//double thresh_v = sqrt(2.0*constants::G*my_planet.get_mass()*(my_parts[0]->get_inverse_radius()-1.0/(my_planet.get_radius()+900e3)));
	cout << "Simulating Particle Transport...\n";

	for (int i=0; i<num_steps; i++)
	{
		if ((i+1) % 10000 == 0)
		{
			cout << i+1 << "\t" << active_parts << endl;
			output_positions("/home/rodney/Documents/coronaTest/data/positions" + to_string(i+1) + ".out");
		}

		for (int j=0; j<num_parts; j++)
		{
			if (my_parts[j]->get_active())
			{
				my_parts[j]->do_timestep(dt, k);
				double r = my_parts[j]->get_radius();
				double temp = common.interpolate(alt_bins, Tn, (r - my_planet.get_radius()));

				if (bg_species.check_collision(r, my_parts[j]->get_total_v(), dt, temp))
				{
					my_parts[j]->do_collision(bg_species.get_collision_target(), bg_species.get_collision_theta());
				}

				// deactivation criteria...need to incorporate into configuration file
				if (my_parts[j]->get_radius() < (my_planet.get_radius() + 900e5) && (my_parts[j]->get_total_v() + v_Obg) < sqrt(2.0*constants::G*my_planet.get_mass()*(my_parts[j]->get_inverse_radius()-1.0/(my_planet.get_radius()+900e5))))
				{
					my_parts[j]->deactivate();
					active_parts--;
				}
			}
		}
	}

	cout << "Number of collisions: " << bg_species.get_num_collisions() << endl;
	cout << "Active particles remaining: " << active_parts << endl;
}
