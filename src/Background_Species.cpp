/*
 * Background_Species.cpp
 *
 *  Created on: Jun 29, 2020
 *      Author: rodney
 */

#include "Background_Species.hpp"

Background_Species::Background_Species() {
	use_temp_profile = false;
	use_dens_profile = false;
	profile_bottom_alt = 0.0;
	profile_top_alt = 0.0;
	num_collisions = 0;
	num_species = 0;
	ref_temp = 0.0;
	ref_height = 0.0;
	ref_g = 0.0;
	collision_target = -1;
	collision_theta = 0.0;
	my_dist = NULL;
}

Background_Species::Background_Species(int n, Planet p, double T, double h, Distribution_MB* dist, Particle* bg_p[], double bg_d[], double bg_s[], string bg_s_f[], string temp_profile_filename, string dens_profile_filename, double profile_bottom, double profile_top) {
	num_collisions = 0;
	num_species = n;
	my_planet = p;
	my_dist = dist;
	ref_temp = T;
	ref_height = h;
	profile_bottom_alt = profile_bottom;
	profile_top_alt = profile_top;
	ref_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+ref_height, 2.0));

	if (temp_profile_filename != "")
	{
		use_temp_profile = true;
	}
	else
	{
		use_temp_profile = false;
	}

	if (dens_profile_filename != "")
	{
		use_dens_profile = true;
	}
	else
	{
		use_dens_profile = false;
	}

	bg_parts.resize(num_species);
	bg_densities.resize(num_species);
	bg_sigma_defaults.resize(num_species);
	bg_sigma_tables.resize(num_species);
	bg_scaleheights.resize(num_species);
	bg_avg_v.resize(num_species);
	for (int i=0; i<num_species; i++)
	{
		bg_parts[i] = bg_p[i];
		bg_densities[i].push_back(bg_d[i]);
		bg_sigma_defaults[i] = bg_s[i];
		bg_sigma_tables[i].resize(2);
		if (bg_s_f[i] != "")
		{
			bg_sigma_defaults[i] = 0.0;
			common::import_csv(bg_s_f[i], bg_sigma_tables[i][0], bg_sigma_tables[i][1]);
		}
		bg_scaleheights[i] = constants::k_b*ref_temp/(bg_parts[i]->get_mass()*ref_g);
		bg_avg_v[i].push_back(sqrt(constants::k_b*ref_temp/bg_parts[i]->get_mass()));
	}
	collision_target = -1;   // set to -1 when no collision happening
	collision_theta = 0.0;

	// read in temperature profile (if available)
	if (use_temp_profile)
	{
		common::import_csv(temp_profile_filename, temp_alt_bins, Tn, Ti, Te);
		int num_alt_bins = temp_alt_bins.size();
		for (int i=0; i<num_species; i++)
		{
			bg_avg_v[i].clear();
			double mass = bg_parts[i]->get_mass();
			for (int j=0; j<num_alt_bins; j++)
			{
				bg_avg_v[i].push_back(sqrt(constants::k_b*Tn[j]/mass));
			}
		}
	}

	// read in density profile (if available)
	if (use_dens_profile)
	{
		for (int i=0; i<num_species; i++)
		{
			bg_densities[i].clear();
		}
		if (num_species == 5)
		{
			common::import_csv(dens_profile_filename, dens_alt_bins, bg_densities[0], bg_densities[1], bg_densities[2], bg_densities[3], bg_densities[4]);
		}
		else if (num_species == 4)
		{
			common::import_csv(dens_profile_filename, dens_alt_bins, bg_densities[0], bg_densities[1], bg_densities[2], bg_densities[3]);
		}
		else if (num_species == 3)
		{
			common::import_csv(dens_profile_filename, dens_alt_bins, bg_densities[0], bg_densities[1], bg_densities[2]);
		}
		else if (num_species == 2)
		{
			common::import_csv(dens_profile_filename, dens_alt_bins, bg_densities[0], bg_densities[1]);
		}
		else if (num_species == 1)
		{
			common::import_csv(dens_profile_filename, dens_alt_bins, bg_densities[0]);
		}
	}

	import_CDF("/home/rodney/Documents/coronaTest/KHARCHENKOCDF.TXT");
}

Background_Species::~Background_Species() {

}

// returns collision energy in eV between particle 1 and particle 2
double Background_Species::calc_collision_e(Particle* p1, Particle* p2)
{
	double energy = 0.0;
	double p1_mass = p1->get_mass();
	double p2_mass = p2->get_mass();
	Matrix<double, 3, 1> p1_v = {p1->get_vx(), p1->get_vy(), p1->get_vz()};
	Matrix<double, 3, 1> p2_v = {p2->get_vx(), p2->get_vy(), p2->get_vz()};
	Matrix<double, 3, 1> vcm;
	vcm = (p1_mass*p1_v.array() + p2_mass*p2_v.array()) / (p1_mass + p2_mass);
	Matrix<double, 3, 1> p1_vcm = p1_v.array() - vcm.array();    // particle 1 c-o-m velocity
	Matrix<double, 3, 1> p2_vcm = p2_v.array() - vcm.array();    // particle 2 c-o-m velocity
	double p1_vcm_tot = sqrt(p1_vcm[0]*p1_vcm[0] + p1_vcm[1]*p1_vcm[1] + p1_vcm[2]*p1_vcm[2]);  // particle 1 c-o-m scalar velocity
	double p2_vcm_tot = sqrt(p2_vcm[0]*p2_vcm[0] + p2_vcm[1]*p2_vcm[1] + p2_vcm[2]*p2_vcm[2]);  // particle 2 c-o-m scalar velocity
	energy = (0.5*p1_mass*p1_vcm_tot*p1_vcm_tot + 0.5*p2_mass*p2_vcm_tot*p2_vcm_tot) / constants::ergev;

	return energy;
}

// calculates new density of background particle based on radial position and scale height
double Background_Species::calc_new_density(double ref_density, double scale_height, double r_moved)
{
	return ref_density*exp(r_moved/scale_height);
}

// check to see if a collision occurred and initialize target particle if so
bool Background_Species::check_collision(Particle* p, double dt)
{
	double r = p->get_radius();
	double my_total_v = p->get_total_v();
	double alt = r - my_planet.get_radius();
	double r_moved = my_planet.get_radius() + ref_height - r;

	// get densities at current location
	vector<double> dens;
	dens.resize(num_species);

	if (use_dens_profile)  // get new density from imported density profile
	{
		for (int i=0; i<num_species; i++)
		{
			dens[i] = get_density(alt, bg_densities[i], bg_parts[i]->get_mass());
		}
	}
	else  // calculate new density based on reference scale height
	{
		for (int i=0; i<num_species; i++)
		{
			dens[i] = calc_new_density(bg_densities[i][0], bg_scaleheights[i], r_moved);
		}
	}

	// look up total cross sections to use for each species, or use default if no table available
	vector<double> total_sig;
	total_sig.resize(num_species);
	for (int i=0; i<num_species; i++)
	{
		double avg_v = 0.0;

		// if default sigma is zero, then table is available, must initialize a particle to get energy
		if (bg_sigma_defaults[i] == 0.0)
		{
			if (use_temp_profile)
			{
				if (alt < profile_bottom_alt)
				{
					avg_v = bg_avg_v[i][0];
				}
				else if (alt > profile_top_alt)
				{
					avg_v = bg_avg_v[i].back();
				}
				else
				{
					avg_v = common::interpolate(temp_alt_bins, bg_avg_v[i], alt);
				}
			}
			else  // use reference temp avg_v
			{
				avg_v = bg_avg_v[i][0];
			}
			my_dist->init_vonly(bg_parts[i], avg_v);

			//double my_mass = p->get_mass();
			//double targ_mass = bg_parts[i]->get_mass();
			//double dm = my_mass*targ_mass / (my_mass + targ_mass);
			//double dv = p->get_total_v() - avg_v;
			//double test_energy = 0.5*dm*dv*dv / constants::ergev;

			// calculate collision energy and look up cross section
			double energy = calc_collision_e(p, bg_parts[i]);
			total_sig[i] = common::interpolate(bg_sigma_tables[i][0], bg_sigma_tables[i][1], energy);

			//cout << "test_energy: " << test_energy << "\t" << "energy: " << energy << "\n";
		}
		else  // just use default sigma if no lookup table available
		{
			total_sig[i] = bg_sigma_defaults[i];
		}
	}

	// determine if test particle collided
	double u = common::get_rand();
	double tau = 0.0;
	for (int i=0; i<num_species; i++)
	{
		tau += 	my_total_v*dt*total_sig[i]*dens[i];
	}
	if (u > exp(-tau))
	{
		num_collisions++;

		// pick target species for collision
		u = common::get_rand();
		double total_dens = 0.0;
		for (int i=0; i<num_species; i++)
		{
			total_dens += dens[i];
		}
		double frac = 0.0;
		collision_target = 0;

		do
		{
			frac += dens[collision_target] / total_dens;
			collision_target++;
		}
		while (u >= frac && collision_target < num_species);

		// subtract the extra added integer, and initialize collision target if necessary
		collision_target--;

		if (bg_sigma_defaults[collision_target] != 0.0)  // particle needs to be initialized
		{
			if (use_temp_profile)
			{
				double avg_v = 0.0;
				if (alt < profile_bottom_alt)
				{
					avg_v = bg_avg_v[collision_target][0];
				}
				else if (alt > profile_top_alt)
				{
					avg_v = bg_avg_v[collision_target].back();
				}
				else
				{
					avg_v = common::interpolate(temp_alt_bins, bg_avg_v[collision_target], alt);
				}
				my_dist->init_vonly(bg_parts[collision_target], avg_v);
			}
			else
			{
				my_dist->init_vonly(bg_parts[collision_target], bg_avg_v[collision_target][0]);
			}
		}
		collision_theta = find_new_theta();
		return true;
	}
	else
	{
		collision_target = -1;
		return false;
	}

	/* experimental
	// determine if test particle collided, and if so, initialize target and pick theta
	double u = common::get_rand();
	double tau = 0.0;
	collision_target = -1;
	int i = 0;
	while (collision_target == -1 && i < num_species)
	{
		// get density at current location
		double dens = 0.0;

		if (use_dens_profile)  // get new density from imported density profile
		{
			dens = get_density(alt, bg_densities[i], bg_parts[i]->get_mass());
		}
		else  // calculate new density based on reference scale height
		{
			dens = calc_new_density(bg_densities[i][0], bg_scaleheights[i], r_moved);
		}

		// look up total cross section to use for this species, or use default if no table available
		double sig = 0.0;
		double avg_v = 0.0;

		// if default sigma is zero, then table is available, must initialize a particle to get energy
		if (bg_sigma_defaults[i] == 0.0)
		{
			if (use_temp_profile)
			{
				if (alt < profile_bottom_alt)
				{
					avg_v = bg_avg_v[i][0];
				}
				else if (alt > profile_top_alt)
				{
					avg_v = bg_avg_v[i].back();
				}
				else
				{
					avg_v = common::interpolate(temp_alt_bins, bg_avg_v[i], alt);
				}
			}
			else  // use reference temp avg_v
			{
				avg_v = bg_avg_v[i][0];
			}
			my_dist->init_vonly(bg_parts[i], avg_v);

			//double my_mass = p->get_mass();
			//double targ_mass = bg_parts[i]->get_mass();
			//double dm = my_mass*targ_mass / (my_mass + targ_mass);
			//double dv = my_total_v - avg_v;
			//double test_energy = 0.5*dm*dv*dv / constants::ergev;
			//double test_energy = 0.5*p->get_mass()*my_total_v*my_total_v / constants::ergev;

			// calculate collision energy and look up cross section
			double energy = calc_collision_e(p, bg_parts[i]);
			sig = common::interpolate(bg_sigma_tables[i][0], bg_sigma_tables[i][1], energy);

			//ofstream outfile;
			//outfile.open("/home/rodney/Documents/coronaTest/energy_compare_HotH_2.out", ios::out | ios::app);
			//outfile << bg_parts[i]->get_name() << "\t" << alt*1e-5 << "km\t" << test_energy << "eV\t" << energy << "eV\n";
			//outfile.close();
		}
		else  // just use default sigma if no lookup table available
		{
			sig = bg_sigma_defaults[i];
		}

		// increase tau and check for collision, move on to next particle if no collision yet
		tau += 	my_total_v*dt*sig*dens;
		if (u > exp(-tau))
		{
			num_collisions++;
			collision_target = i;
			collision_theta = find_new_theta();
			if (bg_sigma_defaults[collision_target] != 0.0)  // particle needs to be initialized
			{
				if (use_temp_profile)
				{
					double avg_v = 0.0;
					if (alt < profile_bottom_alt)
					{
						avg_v = bg_avg_v[collision_target][0];
					}
					else if (alt > profile_top_alt)
					{
						avg_v = bg_avg_v[collision_target].back();
					}
					else
					{
						avg_v = common::interpolate(temp_alt_bins, bg_avg_v[collision_target], alt);
					}
					my_dist->init_vonly(bg_parts[collision_target], avg_v);
				}
				else
				{
					my_dist->init_vonly(bg_parts[collision_target], bg_avg_v[collision_target][0]);
				}
			}
		}
		i++;
	}

	// outside while loop now, check if collision detected
	if (collision_target == -1)
	{
		return false;
	}
	else
	{
		return true;
	}
	*/
}

// scans imported differential scattering CDF for new collision theta
double Background_Species::find_new_theta()
{
	double u = common::get_rand();
	int k = 0;
	while (cdf(k, 0) < u)
	{
		k++;
	}
	return acos(cdf(k, 1));
}

// get density from imported density profile
double Background_Species::get_density(double alt, vector<double> &dens_bins, double targ_mass)
{
	double current_dens = 0.0;

	// if outside of profile boundaries, need to extrapolate using a scale height
	if (alt < profile_bottom_alt)
	{
		double local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+profile_bottom_alt, 2.0));
		current_dens = calc_new_density(dens_bins[0], constants::k_b*common::interpolate(temp_alt_bins, Tn, profile_bottom_alt)/(targ_mass*local_g), profile_bottom_alt - alt);
	}
	else if (alt > profile_top_alt)
	{
		double local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+profile_top_alt, 2.0));
		current_dens = calc_new_density(dens_bins.back(), constants::k_b*common::interpolate(temp_alt_bins, Tn, profile_top_alt)/(targ_mass*local_g), profile_top_alt - alt);
	}
	else
	{
		current_dens = common::interpolate(dens_alt_bins, dens_bins, alt);
	}

	return current_dens;
}

int Background_Species::get_num_collisions()
{
	return num_collisions;
}

Particle* Background_Species::get_collision_target()
{
	return bg_parts[collision_target];
}

double Background_Species::get_collision_theta()
{
	return collision_theta;
}

void Background_Species::import_CDF(string filename)
{
	ifstream infile;
	infile.open(filename);
	string line;

	for (int i=0; i<180; i++)
	{
		getline(infile, line);
		stringstream str(line);

		for (int j=0; j<2; j++)
		{
			str >> cdf(i, j);
		}
	}
	infile.close();
}
