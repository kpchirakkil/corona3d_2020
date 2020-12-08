/*
 * Atmosphere.cpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#include "Atmosphere.hpp"

// construct atmosphere using given parameters
Atmosphere::Atmosphere(int n, int num_to_trace, string trace_output_dir, Planet p, vector<shared_ptr<Particle>> parts, shared_ptr<Distribution> dist, Background_Species bg, int num_EDFs, int EDF_alts[])
{
	num_parts = n;                // number of test particles to track
	num_traced = num_to_trace;    // number of tracked particles to output detailed trace data for
	trace_dir = trace_output_dir;
	active_parts = num_parts;
	my_planet = p;
	my_dist = dist;
	my_parts.resize(num_parts);
	bg_species = bg;

	for (int i=0; i<num_parts; i++)
	{
		my_parts[i] = parts[i];
		my_dist->init(my_parts[i]);
	}

	//initialize stats tracking vectors
	stats_num_EDFs = num_EDFs;
	stats_EDF_alts.resize(stats_num_EDFs);
	stats_EDFs.resize(stats_num_EDFs);
	for (int i=0; i<stats_num_EDFs; i++)
	{
		stats_EDF_alts[i] = EDF_alts[i];
		stats_EDFs[i].resize(1001);
		for (int j=0; j<1001; j++)
		{
			stats_EDFs[i][j] = 0;
		}
	}

	stats_dens_counts.resize(100001);
	for (int i=0; i<100001; i++)
	{
		stats_dens_counts[i] = 0;
	}

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
		string filename = trace_dir + "part" + to_string(traced_parts[i]) + "_collisions.out";
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
			position_file.open(trace_dir + "part" + to_string(traced_parts[i]) + "_positions.out", ios::out | ios::app);
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
void Atmosphere::run_simulation(double dt, int num_steps, double lower_bound, double upper_bound, int print_status_freq, int output_pos_freq, string output_pos_dir)
{
	upper_bound = my_planet.get_radius() + upper_bound;
	lower_bound = my_planet.get_radius() + lower_bound;
	double v_esc_upper = sqrt(2.0 * constants::G * my_planet.get_mass() / upper_bound);
	//double v_esc_lower = sqrt(2.0 * constants::G * my_planet.get_mass() / lower_bound);
	int escape_count = 0;
	int added_particles = 0;  // increment this if re-initializing deactivated particles

	double k = my_planet.get_k_g();
	//double v_Obg = sqrt(8.0*constants::k_b*277.6 / (constants::pi*15.9994*constants::amu));
	cout << "Simulating Particle Transport...\n";

	double vol_at_400 = 2.0*constants::pi/3.0 * (pow(my_planet.get_radius() + 401e5, 3.0) - pow(my_planet.get_radius() + 400e5, 3.0));
    double dens_at_400 = 0.0;
	for (int i=0; i<num_steps; i++)
	{
		if (active_parts == 0)
		{
			break;
		}

		if (print_status_freq > 0 && (i+1) % print_status_freq == 0)
		{
			double hrs = (i+1)*dt/3600.0;
			double min = (hrs - (int)hrs)*60.0;
			double sec = (min - (int)min)*60.0;
			cout << (int)hrs << "h "<< (int)min << "m " << sec << "s " << "\t Active: " << active_parts << "\t Escaped: " << escape_count << "\t Escape fraction: " << (double)escape_count / (double)(num_parts+added_particles) <<endl;

			dens_at_400 = (dt*(1.9097e26/2.0)/(double)(num_parts+added_particles)*stats_dens_counts[400]) / vol_at_400;
			cout << "Average density at 400km: " << dens_at_400 << " per cubic cm\n";
		}

		if (output_pos_freq > 0 && (i+1) % output_pos_freq == 0)
		{
			output_positions(output_pos_dir + "positions" + to_string(i+1) + ".out");
		}

		update_stats();

		if (num_traced > 0)
		{
			output_trace_data();
		}

		for (int j=0; j<num_parts; j++)
		{
			if (my_parts[j]->get_active())
			{
				my_parts[j]->do_timestep(dt, k);

				if (bg_species.check_collision(my_parts[j], dt))
				{
					my_parts[j]->do_collision(bg_species.get_collision_target(), bg_species.get_collision_theta(), i*dt, my_planet.get_radius());
				}

				// deactivation criteria...need to incorporate into configuration file
				//if (my_parts[j]->get_radius() < (my_planet.get_radius() + 900e5) && (my_parts[j]->get_total_v() + v_Obg) < sqrt(2.0*constants::G*my_planet.get_mass()*(my_parts[j]->get_inverse_radius()-1.0/(my_planet.get_radius()+900e5))))
				//{
				//	my_parts[j]->deactivate();
				//	active_parts--;
				//}

				if (my_parts[j]->get_radius() >= upper_bound && my_parts[j]->get_total_v() > v_esc_upper)
				{
					my_parts[j]->deactivate(to_string(i*dt) + "\t\tReached upper bound with escape velocity.\n\n");
					//my_dist->init(my_parts[j]);
					//added_particles++;
					active_parts--;
					escape_count++;
				}
				else if (my_parts[j]->get_radius() <= lower_bound)
				{
					my_parts[j]->deactivate(to_string(i*dt) + "\t\tDropped below lower bound.\n\n");
					//my_dist->init(my_parts[j]);
					//added_particles++;
					active_parts--;
				}
				else if (my_parts[j]->get_total_v() < v_esc_upper)
				{
					my_parts[j]->deactivate(to_string(i*dt) + "\t\tVelocity dropped below upper bound escape velocity.\n\n");
					//my_dist->init(my_parts[j]);
					//added_particles++;
					active_parts--;
				}
			}
		}
	}

	if (num_traced > 0)
	{
		output_collision_data();
	}

	output_stats(dt, 1.9097e26/2, num_parts+added_particles);

	cout << "Number of collisions: " << bg_species.get_num_collisions() << endl;
	cout << "Active particles remaining: " << active_parts << endl;
	cout << "Number of escaped particles: " << escape_count << endl;
	cout << "Fraction of escaped particles: " << (double)escape_count / (double)(num_parts+added_particles) << endl;
}

void Atmosphere::update_stats()
{
	int index = 0;
	double e = 0.0;
	int e_index = 0;

	for (int i=0; i<num_parts; i++)
	{
		if (my_parts[i]->get_active())
		{
			index = (int)(1e-5*(my_parts[i]->get_radius()-my_planet.get_radius()));
			if (index >= 0 && index <= 100000 && my_parts[i]->get_x() > 0.0)
			{
				stats_dens_counts[index]++;
			}
			for (int j=0; j<stats_num_EDFs; j++)
			{
				if (index == stats_EDF_alts[j])
				{
					e = 0.5*my_parts[i]->get_mass()*pow(my_parts[i]->get_total_v(), 2.0)/constants::ergev;
					e_index = (int)(100.0*e);
					if (e_index >= 0 && e_index <= 1000 && my_parts[i]->get_x() > 0.0)
					{
						stats_EDFs[j][e_index]++;
					}
				}
			}
		}
	}
}

void Atmosphere::output_stats(double dt, double rate, int total_parts)
{
	double volume = 0.0;
	double r_in_cm = 0.0;
	double dens = 0.0;
	ofstream dens_out, EDF_out;
	dens_out.open(output_pos_dir + "density1d.out");
	EDF_out.open(output_pos_dir + "EDF_1.out");

	dens_out << "#alt[km]\tdensity[cm-3]\n";
	int size = stats_dens_counts.size();
	for (int i=0; i<size; i++)
	{
		r_in_cm = my_planet.get_radius() + 1e5*(double)i;
		volume = 2.0*constants::pi/3.0 * (pow(r_in_cm+1e5, 3.0) - pow(r_in_cm, 3.0));
	    dens = (dt*rate/(double)total_parts*stats_dens_counts[i]) / volume;

		dens_out << i << "\t\t" << dens << "\n";
	}
	dens_out.close();

	for (int i=0; i<stats_num_EDFs; i++)
	{
		EDF_out.open(output_pos_dir + "EDF_" + to_string(stats_EDF_alts[i]) + "km.out");
		size = stats_EDFs[i].size();
		for (int j=0; j<size; j++)
		{
			EDF_out << (double)i*0.01 << "\t" << stats_EDFs[i][j] << "\n";
		}
		EDF_out.close();
	}
}
