/*
 * Background_Species.cpp
 *
 *  Created on: Jun 29, 2020
 *      Author: rodney
 */

#include "Background_Species.hpp"

Background_Species::Background_Species() {
	use_dens_profile = false;
	num_collisions = 0;
	num_species = 0;
	ref_temp = 0.0;
	ref_height = 0.0;
	ref_g = 0.0;
	collision_target = -1;
	collision_theta = 0.0;
	my_dist = NULL;
}

Background_Species::Background_Species(int n, Planet p, double T, double h, Distribution_MB* dist, Particle* bg_p[], double bg_d[], double bg_s[], string dens_profile_filename) {
	num_collisions = 0;
	num_species = n;
	my_planet = p;
	my_dist = dist;
	ref_temp = T;
	ref_height = h;
	ref_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+ref_height, 2.0));

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
	bg_sigmas.resize(num_species);
	bg_scaleheights.resize(num_species);
	bg_avg_v.resize(num_species);
	for (int i=0; i<num_species; i++)
	{
		bg_parts[i] = bg_p[i];
		bg_densities[i].push_back(bg_d[i]);
		bg_sigmas[i] = bg_s[i];
		bg_scaleheights[i] = constants::k_b*ref_temp/(bg_parts[i]->get_mass()*ref_g);
		bg_avg_v[i] = sqrt(constants::k_b*ref_temp/bg_parts[i]->get_mass());
	}
	collision_target = -1;
	collision_theta = 0.0;

	if (use_dens_profile)
	{
		for (int i=0; i<num_species; i++)
		{
			bg_densities[i].clear();
		}
		if (num_species == 5)
		{
			common.import_csv(dens_profile_filename, alt_bins, bg_densities[0], bg_densities[1], bg_densities[2], bg_densities[3], bg_densities[4]);
		}
		else if (num_species == 4)
		{
			common.import_csv(dens_profile_filename, alt_bins, bg_densities[0], bg_densities[1], bg_densities[2], bg_densities[3]);
		}
		else if (num_species == 3)
		{
			common.import_csv(dens_profile_filename, alt_bins, bg_densities[0], bg_densities[1], bg_densities[2]);
		}
		else if (num_species == 2)
		{
			common.import_csv(dens_profile_filename, alt_bins, bg_densities[0], bg_densities[1]);
		}
		else if (num_species == 1)
		{
			common.import_csv(dens_profile_filename, alt_bins, bg_densities[0]);
		}
	}

	import_CDF("/home/rodney/Documents/coronaTest/KHARCHENKOCDF.TXT");
}

Background_Species::~Background_Species() {

}

// calculates new density of background particle based on radial position and scale height
double Background_Species::calc_new_density(double ref_density, double scale_height, double r_moved)
{
	return ref_density*exp(r_moved/scale_height);
}

// check to see if a collision occurred and set target particle if so
bool Background_Species::check_collision(double r, double v, double dt, double T)
{
	// get densities at current location
	double r_moved = my_planet.get_radius() + ref_height - r;
	std::vector<double> temp_dens;
	temp_dens.resize(num_species);

	if (use_dens_profile)  // get new density from imported density profile
	{
		for (int i=0; i<num_species; i++)
		{
			temp_dens[i] = get_density(r, bg_densities[i]);
		}
	}
	else  // calculate new density based on scale height
	{
		for (int i=0; i<num_species; i++)
		{
			temp_dens[i] = calc_new_density(bg_densities[i][0], bg_scaleheights[i], r_moved);
		}
	}

	// determine if test particle collided
	double u = common.get_rand();
	double tau = 0.0;
	for (int i=0; i<num_species; i++)
	{
		tau += 	v*dt*bg_sigmas[i]*temp_dens[i];
	}
	if (u > exp(-tau))
	{
		num_collisions++;

		// pick target species for collision
		u = common.get_rand();
		double total_dens = 0.0;
		for (int i=0; i<num_species; i++)
		{
			total_dens += temp_dens[i];
		}
		double frac = 0.0;
		collision_target = 0;

		do
		{
			frac += temp_dens[collision_target] / total_dens;
			collision_target++;
		}
		while (u >= frac && collision_target < num_species);

		// subtract the extra added integer, and initialize collision target
		collision_target--;

		if (use_dens_profile)
		{
			double avg_v = sqrt(constants::k_b*T/bg_parts[collision_target]->get_mass());
			my_dist->init_vonly(bg_parts[collision_target], avg_v);
		}
		else
		{
			my_dist->init_vonly(bg_parts[collision_target], bg_avg_v[collision_target]);
		}

		collision_theta = find_new_theta();

		return true;
	}
	else
	{
		collision_target = -1;
		return false;
	}
}

// scans imported differential scattering CDF for new collision theta
double Background_Species::find_new_theta()
{
	double u = common.get_rand();
	int k = 0;
	while (cdf(k, 0) < u)
	{
		k++;
	}
	return acos(cdf(k, 1));
}

// get density from imported density profile
double Background_Species::get_density(double r, vector<double> &dens_bins)
{
	double current_dens = common.interpolate(alt_bins, dens_bins, r - my_planet.get_radius());
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
