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
	stats_EDFs.resize(2);   // index 0 is day side EDFs, 1 is night side
	stats_EDFs[0].resize(stats_num_EDFs);
	stats_EDFs[1].resize(stats_num_EDFs);
	for (int i=0; i<stats_num_EDFs; i++)
	{
		stats_EDF_alts[i] = EDF_alts[i];
		stats_EDFs[0][i].resize(1001);
		stats_EDFs[1][i].resize(1001);
		for (int j=0; j<1001; j++)
		{
			stats_EDFs[0][i][j].resize(1);
			stats_EDFs[1][i][j].resize(1);
			stats_EDFs[0][i][j][0] = 0;
			stats_EDFs[1][i][j][0] = 0;

			/* if using weighting comment out the 4 above lines and uncomment this section
			stats_EDFs[0][i][j].resize(num_parts);
			stats_EDFs[1][i][j].resize(num_parts);
			for (int k=0; k<num_parts; k++)
			{
				stats_EDFs[0][i][j][k] = 0;
				stats_EDFs[1][i][j][k] = 0;
			}
			*/
		}
	}

	stats_dens_counts.resize(2);  // index 0 is day side, 1 is night side
	stats_dens_counts[0].resize(100001);
	stats_dens_counts[1].resize(100001);
	for (int i=0; i<100001; i++)
	{
		stats_dens_counts[0][i].resize(1);
		stats_dens_counts[1][i].resize(1);
		stats_dens_counts[0][i][0] = 0;
		stats_dens_counts[1][i][0] = 0;

		/* if using weighting comment out the 4 above lines and uncomment this section
		stats_dens_counts[0][i].resize(num_parts);
		stats_dens_counts[1][i].resize(num_parts);
		for (int j=0; j<num_parts; j++)
		{
			stats_dens_counts[0][i][j] = 0;
			stats_dens_counts[1][i][j] = 0;
		}
		*/
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
		if (my_parts[i]->is_active())
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
		if (my_parts[i]->is_active())
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
		if (my_parts[traced_parts[i]]->is_active())
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
		if (my_parts[i]->is_active())
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
		if (my_parts[i]->is_active())
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
void Atmosphere::run_simulation(double dt, int num_steps, double lower_bound, double upper_bound, double avg_thermal_v, int print_status_freq, int output_pos_freq, string output_pos_dir, string output_stats_dir)
{
	//int alt_idx = 0;
	//vector<int> collision_counts_by_alt;
	//collision_counts_by_alt.resize(100001);
	//for (int i=0; i<100001; i++)
	//{
	//	collision_counts_by_alt[i] = 0;
	//}

	int thermalized_count = 0;
	int thermal_escape_count = 0;
	int night_escape_count = 0;
	int day_escape_count = 0;
	double v_esc_current = 0.0;
	double wt = 0.0;
	upper_bound = my_planet.get_radius() + upper_bound;
	lower_bound = my_planet.get_radius() + lower_bound;
	double v_esc_upper = sqrt(2.0 * constants::G * my_planet.get_mass() / upper_bound);
	//double v_esc_lower = sqrt(2.0 * constants::G * my_planet.get_mass() / lower_bound);
	int added_parts = 0;  // increment this if re-initializing deactivated particles
	double global_rate = my_dist->get_global_rate();
	double k = my_planet.get_k_g();
	//double v_Obg = sqrt(8.0*constants::k_b*277.6 / (constants::pi*15.9994*constants::amu));

	cout << "Simulating Particle Transport...\n";

	for (int i=0; i<num_steps; i++)
	{
		if (active_parts == 0)
		{
			break;
		}

		// if using weighting, set weight of current particle here
		wt = 1.0;
		//wt = 1.0 / (double)((int)((double)(i+1)*dt) + 1);

		if (print_status_freq > 0 && (i+1) % print_status_freq == 0)
		{
			double hrs = (i+1)*dt/3600.0;
			double min = (hrs - (int)hrs)*60.0;
			double sec = (min - (int)min)*60.0;
			cout << (int)hrs << "h "<< (int)min << "m " << sec << "s " << "\tActive: " << active_parts << "\tDay escape: " << day_escape_count << "\tDay fraction: " << (double)day_escape_count / (double)(num_parts+added_parts) << "\tNight escape: " << night_escape_count <<  "\tNight fraction: " << (double)night_escape_count / (double)(num_parts+added_parts) <<endl;
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
			if (my_parts[j]->is_active())
			{
				my_parts[j]->do_timestep(dt, k);

				if (bg_species.check_collision(my_parts[j], dt))
				{
					//alt_idx = (int)(1e-5*(my_parts[j]->get_radius()-my_planet.get_radius()));
					//if (alt_idx >=0 && alt_idx <=100000)
					//{
					//	collision_counts_by_alt[alt_idx]++;
					//}
					my_parts[j]->do_collision(bg_species.get_collision_target(), bg_species.get_collision_theta(), i*dt, my_planet.get_radius());
				}

				// deactivation criteria from Justin's original Hot O simulation code (must also uncomment v_Obg declaration above to use)
				//if (my_parts[j]->get_radius() < (my_planet.get_radius() + 900e5) && (my_parts[j]->get_total_v() + v_Obg) < sqrt(2.0*constants::G*my_planet.get_mass()*(my_parts[j]->get_inverse_radius()-1.0/(my_planet.get_radius()+900e5))))
				//{
				//	my_parts[j]->deactivate();
				//	active_parts--;
				//}

				// escape velocity at current radius
				v_esc_current = sqrt(2.0 * constants::G * my_planet.get_mass() / my_parts[j]->get_radius());

				if (my_parts[j]->get_radius() >= upper_bound && my_parts[j]->get_total_v() >= v_esc_upper)
				{
					if (my_parts[j]->get_x() > 0.0)
					{
						my_parts[j]->deactivate(to_string(i*dt) + "\t\tReached upper bound on day side with at least escape velocity.\n\n");
						//my_dist->init(my_parts[j]);
						//added_parts++;
						active_parts--;
						day_escape_count++;
					}
					else
					{
						my_parts[j]->deactivate(to_string(i*dt) + "\t\tReached upper bound on night side with at least escape velocity.\n\n");
						//my_dist->init(my_parts[j]);
						//added_parts++;
						active_parts--;
						night_escape_count++;
					}

					if (my_parts[j]->is_thermalized())
					{
						thermal_escape_count++;
					}
				}
				else if (my_parts[j]->get_radius() <= lower_bound)
				{
					my_parts[j]->deactivate(to_string(i*dt) + "\t\tDropped below lower bound.\n\n");
					//my_dist->init(my_parts[j]);
					//added_parts++;
					active_parts--;
				}
				else if (my_parts[j]->get_total_v() < v_esc_current)
				{
					my_parts[j]->deactivate(to_string(i*dt) + "\t\tVelocity dropped below escape velocity.\n\n");
					//my_dist->init(my_parts[j]);
					//added_parts++;
					active_parts--;

					// for tracking thermalized particles instead of deactivating
					if (my_parts[j]->is_thermalized())
					{
						continue;
					}
					else
					{
						my_parts[j]->set_thermalized();
						thermalized_count++;
					}
				}
				/* if using weighting uncomment this section
				else
				{
					my_parts[j]->set_weight(wt);
				}
				*/
			}
		}
	}

	if (num_traced > 0)
	{
		output_collision_data();
	}

	output_stats(dt, (global_rate / 2.0), num_parts+added_parts, output_stats_dir);

	cout << "Number of collisions: " << bg_species.get_num_collisions() << endl;
	cout << "Active particles remaining: " << active_parts << endl;
	cout << "Number of day side escaped particles: " << day_escape_count << endl;
	cout << "Number of night side escaped particles: " << night_escape_count << endl;
	cout << "Total particles spawned: " << num_parts + added_parts << endl;
	cout << "Day side fraction of escaped particles: " << (double)day_escape_count / (double)(num_parts+added_parts) << endl;
	cout << "Night side fraction of escaped particles: " << (double)night_escape_count / (double)(num_parts+added_parts) << endl;
	cout << "Global production rate: " << global_rate << endl;
	cout << "Total thermalized particles: " << thermalized_count << endl;
	cout << "Escaped thermalized particles: " << thermal_escape_count << endl;

	//ofstream outfile;
	//outfile.open(output_stats_dir + "collision_counts.out");
	//outfile << "#alt(km)\tcounts\n";
	//for (int i=0; i<100001; i++)
	//{
	//	outfile << i << "\t\t" << collision_counts_by_alt[i] << "\n";
	//}
	//outfile.close();

	//for (int i=0; i<num_parts; i++)
	//{
	//	cout << my_parts[i]->get_weight() << "\n";
	//}
}

void Atmosphere::update_stats()
{
	int index = 0;
	int prev_index = 0;
	double e = 0.0;
	int e_index = 0;

	for (int i=0; i<num_parts; i++)
	{
		if (my_parts[i]->is_active())
		{
			index = (int)(1e-5*(my_parts[i]->get_radius()-my_planet.get_radius()));
			prev_index = (int)(1e-5*(my_parts[i]->get_previous_radius()-my_planet.get_radius()));
			if (index >= 0 && index <= 100000)// && prev_index == (index-1))
			{
				if (my_parts[i]->get_x() > 0.0)
				{
					// if using weighting use second line below and comment out other
					stats_dens_counts[0][index][0]++;
					//stats_dens_counts[0][index][i]++;
				}
				else
				{
					// if using weighting use second line below and comment out other
					stats_dens_counts[1][index][0]++;
					//stats_dens_counts[1][index][i]++;
				}
			}

			for (int j=0; j<stats_num_EDFs; j++)
			{
				if (index == stats_EDF_alts[j])
				{
					e = my_parts[i]->get_energy_in_eV();
					e_index = (int)(100.0*e);
					if (e_index >= 0 && e_index <= 1000)// && prev_index == (index-1))
					{
						if (my_parts[i]->get_x() > 0.0)
						{
							// if using weighting use second line below and comment out other
							stats_EDFs[0][j][e_index][0]++;
							//stats_EDFs[0][j][e_index][i]++;
						}
						else
						{
							// if using weighting use second line below and comment out other
							stats_EDFs[1][j][e_index][0]++;
							//stats_EDFs[1][j][e_index][i]++;
						}
					}
				}
			}
		}
	}
}

void Atmosphere::output_stats(double dt, double rate, int total_parts, string output_dir)
{
	vector<double> normed_EDF_day;
	vector<double> normed_EDF_night;
	double volume = 0.0;
	double r_in_cm = 0.0;
	double dens_day = 0.0;
	double dens_night = 0.0;
	double weighted_sum_day = 0.0;
	double weighted_sum_night = 0.0;
	ofstream dens_day_out, dens_night_out, EDF_day_out, EDF_night_out;
	dens_day_out.open(output_dir + "density1d_day.out");
	dens_night_out.open(output_dir + "density1d_night.out");

	dens_day_out << "#alt[km]\tdensity[cm-3]\n";
	dens_night_out << "#alt[km]\tdensity[cm-3]\n";
	int size = stats_dens_counts[0].size();
	for (int i=0; i<size; i++)
	{
		weighted_sum_day = 0.0;
		weighted_sum_night = 0.0;
		r_in_cm = my_planet.get_radius() + 1e5*(double)i;
		volume = 2.0*constants::pi/3.0 * (pow(r_in_cm+1e5, 3.0) - pow(r_in_cm, 3.0));

		weighted_sum_day = (double)stats_dens_counts[0][i][0];
		weighted_sum_night = (double)stats_dens_counts[1][i][0];

		/* if using weighting, uncomment this section and comment out above 2 lines
		for (int j=0; j<num_parts; j++)
		{
			weighted_sum_day += (double)stats_dens_counts[0][i][j] * my_parts[j]->get_weight();
			weighted_sum_night += (double)stats_dens_counts[1][i][j] * my_parts[j]->get_weight();
		}
		*/
	    dens_day = (dt*rate/(double)total_parts*weighted_sum_day) / volume;
	    dens_night = (dt*rate/(double)total_parts*weighted_sum_night) / volume;

		dens_day_out << i << "\t\t" << dens_day << "\n";
		dens_night_out << i << "\t\t" << dens_night << "\n";
	}
	dens_day_out.close();
	dens_night_out.close();

	for (int i=0; i<stats_num_EDFs; i++)
	{
		EDF_day_out.open(output_dir + "EDF_day_" + to_string(stats_EDF_alts[i]) + "km.out");
		EDF_night_out.open(output_dir + "EDF_night_" + to_string(stats_EDF_alts[i]) + "km.out");
		EDF_day_out << "#energy[eV]\tdistribution[cm-3 eV-1]\n";
		EDF_night_out << "#energy[eV]\tdistribution[cm-3 eV-1]\n";
		weighted_sum_day = 0.0;
		weighted_sum_night = 0.0;

		size = stats_EDFs[0][i].size();
		normed_EDF_day.clear();
		normed_EDF_night.clear();
		normed_EDF_day.reserve(size);
		normed_EDF_night.reserve(size);
		for (int j=0; j<size; j++)
		{
			normed_EDF_day[j] = 0.0;
			normed_EDF_night[j] = 0.0;

			normed_EDF_day[j] = (double)stats_EDFs[0][i][j][0];
			normed_EDF_night[j] = (double)stats_EDFs[1][i][j][0];
			/* if using weighting, uncomment this section and comment out above 2 lines
			for (int k=0; k<num_parts; k++)
			{
				normed_EDF_day[j] += (double)stats_EDFs[0][i][j][k] * my_parts[k]->get_weight();
				normed_EDF_night[j]  += (double)stats_EDFs[1][i][j][k] * my_parts[k]->get_weight();
			}
			*/
			weighted_sum_day += normed_EDF_day[j];
			weighted_sum_night += normed_EDF_night[j];
		}
		for (int j=0; j<size; j++)
		{
			normed_EDF_day[j] = normed_EDF_day[j] / weighted_sum_day;
			normed_EDF_night[j] = normed_EDF_night[j] / weighted_sum_night;
		}

		r_in_cm = 1e5*(double)stats_EDF_alts[i] + my_planet.get_radius();
		volume = 2.0*constants::pi/3.0 * (pow(r_in_cm+1e5, 3.0) - pow(r_in_cm, 3.0));
		dens_day = (dt*rate/(double)total_parts*weighted_sum_day) / volume;
		dens_night = (dt*rate/(double)total_parts*weighted_sum_night) / volume;

		for (int j=0; j<size; j++)
		{
			EDF_day_out << (double)j*0.01 << "\t\t" << normed_EDF_day[j]*(dens_day/0.01) << "\n";
			EDF_night_out << (double)j*0.01 << "\t\t" << normed_EDF_night[j]*(dens_night/0.01) << "\n";
		}
		EDF_day_out.close();
		EDF_night_out.close();
	}
}
