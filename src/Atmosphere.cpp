/*
 * Atmosphere.cpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#include "Atmosphere.hpp"

// construct atmosphere using given parameters
Atmosphere::Atmosphere(int n, int num_to_trace, Planet p, vector<Particle*> parts, Distribution* dist, Background_Species bg)
{
	num_parts = n;                // number of test particles to track
	num_traced = num_to_trace;    // number of tracked particles to output detailed trace data for
	active_parts = num_parts;
	my_planet = p;
	my_dist = dist;
	my_parts.resize(num_parts);
	bg_species = bg;
	shell_active = false;
	shell_bottom = 0.0;
	shell_top = 0.0;
	shell_enter_top = 0;
	shell_enter_bottom = 0;
	shell_exit_top = 0;
	shell_exit_bottom = 0;
	shell_numvelbins = 0;
	shell_velbinwidth = 0.0;

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

	// pick trace particles if any
	if (num_traced > 0)
	{
		traced_parts.resize(num_traced);
		for (int i=0; i<num_traced; i++)
		{
			traced_parts[i] = common::get_rand_int(0, num_parts-1);
			my_parts[traced_parts[i]]->set_traced();
		}
	}
}

Atmosphere::~Atmosphere() {

}

// initialize a shell of atmosphere between bottom_r and top_r for data collection
void Atmosphere::init_shell(double bottom_r, double top_r, int num_bins, double bin_width)
{
	shell_active = true;
	shell_bottom = bottom_r;
	shell_top = top_r;
	shell_enter_top = 0;
	shell_enter_bottom = 0;
	shell_exit_top = 0;
	shell_exit_bottom = 0;
	shell_numvelbins = num_bins;
	shell_velbinwidth = bin_width;
	shell_velbins.clear();
	shell_velbins.resize(shell_numvelbins);
}

void Atmosphere::output_shell_data(string folder_path)
{
	ofstream vel_outfile, info_outfile;
	vel_outfile.open(folder_path + "shell_veldistro.out");
	info_outfile.open(folder_path + "shell_info.out");

	vel_outfile << shell_velbinwidth << "\n";
	vel_outfile << shell_numvelbins << "\n";

	for (int i=0; i<shell_numvelbins; i++)
	{
		vel_outfile << shell_velbins[i] << '\n';
	}
	vel_outfile.close();

	info_outfile << "Shell bottom: \t\t" << to_string(shell_bottom) << "\n";
	info_outfile << "Shell top: \t\t" << to_string(shell_top) << "\n";
	info_outfile << "Enter bottom: \t\t" << to_string(shell_enter_bottom) << "\n";
	info_outfile << "Enter top: \t\t" << to_string(shell_enter_top) << "\n";
	info_outfile << "Exit bottom: \t\t" << to_string(shell_exit_bottom) << "\n";
	info_outfile << "Exit top: \t\t" << to_string(shell_exit_top) << "\n";
	info_outfile.close();

	shell_active = false;
}

// writes single-column output file of altitude bin counts using active particles
// bin_width is in cm; first 2 numbers in output file are bin_width and num_bins
void Atmosphere::output_altitude_distro(double bin_width, string datapath)
{
	ofstream outfile;
	outfile.open(datapath);
	double alt = 0.0;           // altitude of particle [cm]
	int nb = 0;                 // bin number

	double max_radius = 0.0;
	for (int i=0; i<num_parts; i++)
	{
		if (my_parts[i]->get_active())
		{
			if (my_parts[i]->get_radius() > max_radius)
			{
				max_radius = my_parts[i]->get_radius();
			}
		}
	}
	int num_bins = (int)((max_radius - my_planet.get_radius()) / bin_width) + 10;
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

	outfile << bin_width << "\n";
	outfile << num_bins << "\n";

	for (int i=0; i<num_bins; i++)
	{
		outfile << abins[i] << "\n";
	}
	outfile.close();
}

void Atmosphere::output_collision_data()
{
	for (int i=0; i<num_traced; i++)
	{
		string filename = "/home/rodney/Documents/coronaTest/trace_data/part" + to_string(traced_parts[i]) + "_collisions.out";
		my_parts[traced_parts[i]]->dump_collision_log(filename);
	}
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

// output test particle trace data for selected particles
void Atmosphere::output_trace_data()
{
	for (int i=0; i<num_traced; i++)
	{
		if (my_parts[traced_parts[i]]->get_active())
		{
			ofstream position_file;
			position_file.open("/home/rodney/Documents/coronaTest/trace_data/part" + to_string(traced_parts[i]) + "_positions.out", ios::out | ios::app);
			position_file << setprecision(10) << my_parts[traced_parts[i]]->get_x() << '\t';
			position_file << setprecision(10) << my_parts[traced_parts[i]]->get_y() << '\t';
			position_file << setprecision(10) << my_parts[traced_parts[i]]->get_z() << '\n';
			position_file.close();
		}
	}
}

// writes single-column output file of velocity bin counts using active particles
// bin_width is in cm/s; first 2 numbers in output file are bin_width and num_bins
void Atmosphere::output_velocity_distro(double bin_width, string datapath)
{
	ofstream outfile;
	outfile.open(datapath);
	double v = 0.0;             // velocity magnitude [cm/s]
	int nb = 0;                 // bin number

	double max_v = 0.0;
	for (int i=0; i<num_parts; i++)
	{
		if (my_parts[i]->get_active())
		{
			double total_v = my_parts[i]->get_total_v();
			if (total_v > max_v)
			{
				max_v = total_v;
			}
		}
	}
	int num_bins = (int)((max_v / bin_width) + 10);
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

	outfile << bin_width << "\n";
	outfile << num_bins << "\n";

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
	double v_Obg = sqrt(8.0*constants::k_b*277.6 / (constants::pi*15.9994*constants::amu));
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

		if (num_traced > 0)
		{
			output_trace_data();
		}

		for (int j=0; j<num_parts; j++)
		{
			if (my_parts[j]->get_active())
			{
				my_parts[j]->do_timestep(dt, k);

				if (bg_species.check_collision(my_parts[j]->get_radius(), my_parts[j]->get_total_v(), dt))
				{
					my_parts[j]->do_collision(bg_species.get_collision_target(), bg_species.get_collision_theta(), i*dt, my_planet.get_radius());
				}

				// deactivation criteria...need to incorporate into configuration file
				if (my_parts[j]->get_radius() < (my_planet.get_radius() + 900e5) && (my_parts[j]->get_total_v() + v_Obg) < sqrt(2.0*constants::G*my_planet.get_mass()*(my_parts[j]->get_inverse_radius()-1.0/(my_planet.get_radius()+900e5))))
				{
					my_parts[j]->deactivate();
					active_parts--;
				}
			}
		}
		// update shell tracking data if necessary
		if (shell_active)
		{
			update_shell_data();
		}
	}

	if (num_traced > 0)
	{
		output_collision_data();
	}

	cout << "Number of collisions: " << bg_species.get_num_collisions() << endl;
	cout << "Active particles remaining: " << active_parts << endl;
}

// update shell numbers (velocity distro, etc.)
void Atmosphere::update_shell_data()
{
	int bin_num = 0;
	double v = 0.0;

	for (int i=0; i<num_parts; i++)
	{
		if (my_parts[i]->get_active())
		{
			double r = my_parts[i]->get_radius();
			double prev_r = my_parts[i]->get_previous_radius();
			if (r > shell_bottom && r < shell_top)
			{
				if (prev_r < shell_bottom)
				{
					shell_enter_bottom++;
					v = my_parts[i]->get_total_v();
					bin_num = (int)(v / shell_velbinwidth);
					shell_velbins[bin_num]++;
				}
				else if (prev_r > shell_top)
				{
					shell_enter_top++;
					v = my_parts[i]->get_total_v();
					bin_num = (int)(v / shell_velbinwidth);
					shell_velbins[bin_num]++;
				}
			}
			else if (prev_r > shell_bottom && prev_r < shell_top)
			{
				if (r < shell_bottom)
				{
					shell_exit_bottom++;
				}
				else if (r > shell_top)
				{
					shell_exit_top++;
				}
			}
		}
	}
}
