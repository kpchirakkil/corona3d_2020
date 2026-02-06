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
	num_inelastic_collisions = 0;
	num_species = 0;
	ref_temp = 0.0;
	ref_height = 0.0;
	ref_g = 0.0;
	collision_target = -1;
	collision_theta = 0.0;
	my_dist = NULL;
	last_outcome = {false, false, -1, 0.0, 0.0};
}

Background_Species::Background_Species(int num_parts, string config_files[], Planet p, double ref_T, double ref_h, string temp_profile_filename, string dens_profile_filename, double profile_bottom, double profile_top)
{
	num_collisions = 0;
	num_inelastic_collisions = 0;
	last_outcome = {false, false, -1, 0.0, 0.0};
	num_species = num_parts;
	my_planet = p;
	ref_temp = ref_T;
	ref_height = ref_h;
	my_dist = make_shared<Distribution_MB>(my_planet, ref_height, ref_temp);
	profile_bottom_alt = profile_bottom;
	profile_top_alt = profile_top;
	collision_target = -1;   // set to -1 when no collision happening
	collision_theta = 0.0;
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
	dens_interp.resize(num_species);
	bg_sigma_defaults.resize(num_species);
	bg_sigma_tables.resize(num_species);
	sigma_interp.resize(num_species);
	bg_scaleheights.resize(num_species);
	bg_avg_v.resize(num_species);
	avg_v_interp.resize(num_species);
	diff_sigma_energies.resize(num_species);
	diff_sigma_CDFs.resize(num_species);
	enable_inelastic.resize(num_species, false);
	sigma_total_interp.resize(num_species);
	elastic_frac_interp.resize(num_species);
	avg_eloss_interp.resize(num_species);
	inelastic_CDFs.resize(num_species);
	for (int i=0; i<num_species; i++)
	{
		int num_energies = 0;
		int energies_index = 0;
		bg_sigma_tables[i].resize(2);
		diff_sigma_CDFs[i].resize(2);

		ifstream infile;
		infile.open(config_files[i]);
		if (!infile.good())
		{
			cout << "Background species configuration file " + to_string(i+1) + " not found!\n";
			exit(1);
		}
		string line, param, val;
		vector<string> parameters;
		vector<string> values;
		int num_params = 0;

				while (getline(infile, line))
                {
                        if (line.empty() || line[0] == '#' || std::all_of(line.begin(), line.end(), ::isspace))
                        {
                                continue;
                        }
			else
			{
				stringstream str(line);
				str >> param >> val;
				parameters.push_back(param);
				values.push_back(val);
				num_params++;
				param = "";
				val = "";
			}
		}
		infile.close();

		for (int j=0; j<num_params; j++)
		{
			if (parameters[j] == "type")
			{
				bg_parts[i] = set_particle_type(values[j]);
			}
			else if (parameters[j] == "ref_dens")
			{
				bg_densities[i].push_back(stod(values[j]));
			}
			else if (parameters[j] == "total_sigma_default")
			{
				bg_sigma_defaults[i] = stod(values[j]);
			}
			else if (parameters[j] == "total_sigma_file")
			{
				if (values[j] != "")
				{
					bg_sigma_defaults[i] = 0.0;
					common::import_csv(values[j], bg_sigma_tables[i][0], bg_sigma_tables[i][1]);
					sigma_interp[i] = make_shared<Interpolator>(bg_sigma_tables[i][0], bg_sigma_tables[i][1]);
				}
			}
			else if (parameters[j] == "num_diff_energies")
			{
				num_energies = stoi(values[j]);
				energies_index = j+1;
				diff_sigma_CDFs[i].resize(num_energies);
			}
			else if (parameters[j] == "enable_inelastic")
			{
				enable_inelastic[i] = (values[j] == "true");
			}
			else if (parameters[j] == "total_sigma_file_total")
			{
				if (values[j] != "")
				{
					vector<vector<double>> table(2);
					common::import_csv(values[j], table[0], table[1]);
					sigma_total_interp[i] = make_shared<Interpolator>(table[0], table[1]);
				}
			}
			else if (parameters[j] == "elastic_fraction_file")
			{
				if (values[j] != "")
				{
					vector<vector<double>> table(2);
					common::import_csv(values[j], table[0], table[1]);
					elastic_frac_interp[i] = make_shared<Interpolator>(table[0], table[1]);
				}
			}
			else if (parameters[j] == "avg_energy_loss_file")
			{
				if (values[j] != "")
				{
					vector<vector<double>> table(2);
					common::import_csv(values[j], table[0], table[1]);
					avg_eloss_interp[i] = make_shared<Interpolator>(table[0], table[1]);
				}
			}
		}

		// Load inelastic DCS files if inelastic is enabled
		if (enable_inelastic[i])
		{
			inelastic_CDFs[i].resize(num_energies);
			for (int j=0; j<num_params; j++)
			{
				// Match energyN_inelastic_file pattern
				string p = parameters[j];
				if (p.length() > 15 && p.substr(p.length() - 15) == "_inelastic_file")
				{
					// Extract energy index from parameter name (e.g., "energy1_inelastic_file" -> 0)
					string num_str = p.substr(6, p.length() - 21);
					int eidx = stoi(num_str) - 1;  // 0-based

					if (eidx >= 0 && eidx < num_energies)
					{
						inelastic_CDFs[i][eidx].resize(2);
						vector<vector<double>> inel_PDF(2);
						common::import_csv(values[j], inel_PDF[0], inel_PDF[1]);
						make_new_inelastic_CDF(i, eidx, inel_PDF[0], inel_PDF[1]);
					}
				}
			}
			// Verify all inelastic CDFs were loaded
			int loaded_count = 0;
			for (int k=0; k<num_energies; k++)
			{
				if (inelastic_CDFs[i][k].size() >= 2 && !inelastic_CDFs[i][k][0].empty())
					loaded_count++;
			}
			if (loaded_count < num_energies)
			{
				cout << "ERROR: enable_inelastic=true for species " << bg_parts[i]->get_name()
				     << " but only " << loaded_count << " of " << num_energies
				     << " inelastic DCS files were loaded. Check config for missing energyN_inelastic_file entries." << endl;
				exit(1);
			}
			cout << "Inelastic collisions enabled for species " << bg_parts[i]->get_name()
			     << " (" << loaded_count << "/" << num_energies << " inelastic CDFs loaded)" << endl;
		}

		bg_scaleheights[i].push_back(constants::k_b*ref_temp/(bg_parts[i]->get_mass()*ref_g));
		bg_avg_v[i].push_back(sqrt(constants::k_b*ref_temp/bg_parts[i]->get_mass()));

		for (int j=0; j<num_energies; j++)
		{
			diff_sigma_CDFs[i][j].resize(2);
			vector<vector<double>> diff_sigma_PDF;
			diff_sigma_PDF.resize(2);
			diff_sigma_energies[i].push_back(stod(values[energies_index + j]));
			common::import_csv(values[energies_index + num_energies + j], diff_sigma_PDF[0], diff_sigma_PDF[1]);
			make_new_CDF(i, j, diff_sigma_PDF[0], diff_sigma_PDF[1]);
		}
	}

	// read in temperature profile (if available) and set avg_v for each alt bin for each species
	if (use_temp_profile)
	{
		common::import_csv(temp_profile_filename, temp_alt_bins, Tn, Ti, Te);
		Tn_interp = make_shared<Interpolator>(temp_alt_bins, Tn);
		Ti_interp = make_shared<Interpolator>(temp_alt_bins, Ti);
		Te_interp = make_shared<Interpolator>(temp_alt_bins, Te);
		int num_alt_bins = temp_alt_bins.size();
		for (int i=0; i<num_species; i++)
		{
			// clear avg_v based on ref_temp and populate with values based on temp profile
			bg_avg_v[i].clear();
			double mass = bg_parts[i]->get_mass();
			for (int j=0; j<num_alt_bins; j++)
			{
				bg_avg_v[i].push_back(sqrt(constants::k_b*Tn[j]/mass));
			}

			// generate interpolator for bg_avg_v for each species
			avg_v_interp[i] = make_shared<Interpolator>(temp_alt_bins, bg_avg_v[i]);
		}
	}

	// read in density profile (if available)
	if (use_dens_profile)
	{
		for (int i=0; i<num_species; i++)
		{
			// clear out default densities and scale heights in order to use profile derived values
			bg_densities[i].clear();
			bg_scaleheights[i].clear();
			bg_scaleheights[i].resize(2);
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

		double bottom_local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+profile_bottom_alt, 2.0));
		double top_local_g = (constants::G * my_planet.get_mass()) / (pow(my_planet.get_radius()+profile_top_alt, 2.0));
		for (int i=0; i<num_species; i++)
		{
			// generate interpolator for each density profile
			dens_interp[i] = make_shared<Interpolator>(dens_alt_bins, bg_densities[i]);

			// calc top and bottom scale height to be used for extrapolating densities
			bg_scaleheights[i][0] = constants::k_b*Tn_interp->loglinterp(profile_bottom_alt)/(bg_parts[i]->get_mass()*bottom_local_g);
			bg_scaleheights[i][1] = constants::k_b*Tn_interp->loglinterp(profile_top_alt)/(bg_parts[i]->get_mass()*top_local_g);
		}
	}
}


Background_Species::~Background_Species() {

}

// returns collision energy in eV between particle 1 and particle 2
double Background_Species::calc_collision_e(shared_ptr<Particle> p1, shared_ptr<Particle> p2)
{
	double e = 0.0;
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
	e = (0.5*p1_mass*p1_vcm_tot*p1_vcm_tot + 0.5*p2_mass*p2_vcm_tot*p2_vcm_tot) / constants::ergev;

	return e;
}

// calculates new density of background particle based on radial position and scale height
double Background_Species::calc_new_density(double ref_density, double scale_height, double r_moved)
{
	return ref_density*exp(r_moved/scale_height);
}

// check to see if a collision occurred and initialize target particle if so
bool Background_Species::check_collision(shared_ptr<Particle> p, double dt)
{
	vector<double> energy;
	energy.resize(num_species);
	double r = p->get_radius();
	double my_total_v = p->get_total_v();
	double alt = r - my_planet.get_radius();
	double r_moved = my_planet.get_radius() + ref_height - r;

	// OPTIMIZATION: Skip collision check at high altitudes (>1000 km)
	// where atmospheric density is negligible and collisions essentially never occur.
	// This significantly speeds up particle tracking in the exosphere.
	if (alt > 1000e5) {  // 1000 km in cm
		collision_target = -1;
		return false;
	}

	// get densities at current location
	vector<double> dens;
	dens.resize(num_species);

	if (use_dens_profile)  // get new density from imported density profile
	{
		for (int i=0; i<num_species; i++)
		{
			dens[i] = get_density(alt, i);
		}
	}
	else  // calculate new density based on reference scale height
	{
		for (int i=0; i<num_species; i++)
		{
			dens[i] = calc_new_density(bg_densities[i][0], bg_scaleheights[i][0], r_moved);
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
					avg_v = avg_v_interp[i]->loglinterp(alt);
				}
			}
			else  // use reference temp avg_v
			{
				avg_v = bg_avg_v[i][0];
			}
			my_dist->init_vonly(bg_parts[i], avg_v);

			// calculate collision energy and look up cross section
			energy[i] = calc_collision_e(p, bg_parts[i]);
			if (enable_inelastic[i] && sigma_total_interp[i])
				total_sig[i] = sigma_total_interp[i]->linterp(energy[i]);
			else
				total_sig[i] = sigma_interp[i]->linterp(energy[i]);
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

		// pick target species for collision, weighted by dens[i] * sigma[i]
		u = common::get_rand();
		double total_weight = 0.0;
		for (int i=0; i<num_species; i++)
		{
			total_weight += dens[i] * total_sig[i];
		}
		double frac = 0.0;
		collision_target = 0;

		do
		{
			frac += (dens[collision_target] * total_sig[collision_target]) / total_weight;
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
					avg_v = avg_v_interp[collision_target]->loglinterp(alt);
				}
				my_dist->init_vonly(bg_parts[collision_target], avg_v);
			}
			else
			{
				my_dist->init_vonly(bg_parts[collision_target], bg_avg_v[collision_target][0]);
			}
			energy[collision_target] = calc_collision_e(p, bg_parts[collision_target]);
		}
		last_outcome.occurred = true;
		last_outcome.target_index = collision_target;
		last_outcome.is_inelastic = false;
		last_outcome.delta_E_eV = 0.0;

		if (enable_inelastic[collision_target] && elastic_frac_interp[collision_target])
		{
			double f_el = elastic_frac_interp[collision_target]->linterp(energy[collision_target]);
			f_el = max(0.0, min(1.0, f_el));

			if (common::get_rand() > f_el)
			{
				// INELASTIC
				last_outcome.is_inelastic = true;
				last_outcome.delta_E_eV = avg_eloss_interp[collision_target]->linterp(energy[collision_target]);
				last_outcome.theta = find_new_theta_inelastic(collision_target, energy[collision_target]);
				num_inelastic_collisions++;
			}
			else
			{
				// ELASTIC
				last_outcome.theta = find_new_theta(collision_target, energy[collision_target]);
			}
		}
		else
		{
			// No inelastic data — elastic only (unchanged behavior)
			last_outcome.theta = find_new_theta(collision_target, energy[collision_target]);
		}

		collision_theta = last_outcome.theta;  // keep backwards-compatible member
		return true;
	}
	else
	{
		collision_target = -1;
		return false;
	}
}

// scans imported differential scattering CDF for new collision theta
double Background_Species::find_new_theta(int part_index, double energy)
{
	// get energy index
	int energy_index = 0;
	int num_energies = diff_sigma_energies[part_index].size();
	if (energy <= diff_sigma_energies[part_index][0])
	{
		energy_index = 0;
	}
	else if (energy >= diff_sigma_energies[part_index].back())
	{
		energy_index = num_energies - 1;
	}
	else
	{
		double difference = INFINITY;
		for (int i=0; i<num_energies; i++)
		{
			double new_diff = abs(energy - diff_sigma_energies[part_index][i]);
			if (new_diff < difference)
			{
				difference = new_diff;
				energy_index = i;
			}
		}
	}

		// search CDF for angle
        double u = common::get_rand();
        // clamp probability to [0,1] in case of rounding errors
        u = std::max(0.0, std::min(u, 1.0));
        int k = 0;
        while ((k + 1 < (int)diff_sigma_CDFs[part_index][energy_index][0].size()) &&
               diff_sigma_CDFs[part_index][energy_index][0][k] < u)
        {
                k++;
        }
        return diff_sigma_CDFs[part_index][energy_index][1][k];
}

// scans imported inelastic differential scattering CDF for new collision theta
double Background_Species::find_new_theta_inelastic(int part_index, double energy)
{
	// get energy index (uses same energy grid as elastic DCS)
	int energy_index = 0;
	int num_energies = diff_sigma_energies[part_index].size();
	if (energy <= diff_sigma_energies[part_index][0])
	{
		energy_index = 0;
	}
	else if (energy >= diff_sigma_energies[part_index].back())
	{
		energy_index = num_energies - 1;
	}
	else
	{
		double difference = INFINITY;
		for (int i=0; i<num_energies; i++)
		{
			double new_diff = abs(energy - diff_sigma_energies[part_index][i]);
			if (new_diff < difference)
			{
				difference = new_diff;
				energy_index = i;
			}
		}
	}

	// Guard: fall back to elastic DCS if inelastic CDF data is missing
	if (energy_index >= (int)inelastic_CDFs[part_index].size() ||
	    inelastic_CDFs[part_index][energy_index].size() < 2 ||
	    inelastic_CDFs[part_index][energy_index][0].empty())
	{
		return find_new_theta(part_index, energy);
	}

	// search inelastic CDF for angle
	double u = common::get_rand();
	u = std::max(0.0, std::min(u, 1.0));
	int k = 0;
	while ((k + 1 < (int)inelastic_CDFs[part_index][energy_index][0].size()) &&
	       inelastic_CDFs[part_index][energy_index][0][k] < u)
	{
		k++;
	}
	return inelastic_CDFs[part_index][energy_index][1][k];
}

// get density from imported density profile
double Background_Species::get_density(double alt, int index)
{
	double current_dens = 0.0;

	// if outside of profile boundaries, need to extrapolate using a scale height
	if (alt < profile_bottom_alt)
	{
		current_dens = calc_new_density(bg_densities[index][0], bg_scaleheights[index][0], profile_bottom_alt - alt);
	}
	else if (alt > profile_top_alt)
	{
		current_dens = calc_new_density(bg_densities[index].back(), bg_scaleheights[index][1], profile_top_alt - alt);
	}
	else
	{
		current_dens = dens_interp[index]->loglinterp(alt);
	}

	return current_dens;
}

int Background_Species::get_num_collisions()
{
	return num_collisions;
}

shared_ptr<Particle> Background_Species::get_collision_target()
{
	return bg_parts[collision_target];
}

double Background_Species::get_collision_theta()
{
	return collision_theta;
}

CollisionOutcome Background_Species::get_last_outcome()
{
	return last_outcome;
}

int Background_Species::get_num_inelastic_collisions()
{
	return num_inelastic_collisions;
}

// make a new differential cross section CDF and store at diff_sigma_CDFs[part_index][energy_index]
void Background_Species::make_new_CDF(int part_index, int energy_index, vector<double> &angle, vector<double> &sigma)
{
	int num_angles = angle.size();
	diff_sigma_CDFs[part_index][energy_index][0].resize(num_angles);
	diff_sigma_CDFs[part_index][energy_index][1].resize(num_angles);

	double sig_total = 0.0;
	for (int i=0; i<num_angles; i++)
	{
		diff_sigma_CDFs[part_index][energy_index][1][i] = angle[i] * (constants::pi / 180.0);
		sigma[i] = sigma[i] * sin(angle[i]*constants::pi/180.0);
		sig_total = sig_total + sigma[i];
	}
	for (int i=0; i<num_angles; i++)
	{
		if (i == 0)
		{
			diff_sigma_CDFs[part_index][energy_index][0][i] = sigma[i] / sig_total;
		}
		else
		{
			diff_sigma_CDFs[part_index][energy_index][0][i] = (sigma[i] / sig_total) + diff_sigma_CDFs[part_index][energy_index][0][i-1];
		}
	}
}

// make a new inelastic differential cross section CDF and store at inelastic_CDFs[part_index][energy_index]
void Background_Species::make_new_inelastic_CDF(int part_index, int energy_index, vector<double> &angle, vector<double> &sigma)
{
	int num_angles = angle.size();
	inelastic_CDFs[part_index][energy_index][0].resize(num_angles);
	inelastic_CDFs[part_index][energy_index][1].resize(num_angles);

	double sig_total = 0.0;
	for (int i=0; i<num_angles; i++)
	{
		inelastic_CDFs[part_index][energy_index][1][i] = angle[i] * (constants::pi / 180.0);
		sigma[i] = sigma[i] * sin(angle[i]*constants::pi/180.0);
		sig_total = sig_total + sigma[i];
	}
	for (int i=0; i<num_angles; i++)
	{
		if (i == 0)
		{
			inelastic_CDFs[part_index][energy_index][0][i] = sigma[i] / sig_total;
		}
		else
		{
			inelastic_CDFs[part_index][energy_index][0][i] = (sigma[i] / sig_total) + inelastic_CDFs[part_index][energy_index][0][i-1];
		}
	}
}

//subroutine to set particle types
shared_ptr<Particle> Background_Species::set_particle_type(string type)
{
	shared_ptr<Particle> p;

	if (type == "H")
	{
		p = make_shared<Particle_H>();
	}
	else if (type == "O")
	{
		p = make_shared<Particle_O>();
	}
	else if (type == "N2")
	{
		p = make_shared<Particle_N2>();
	}
	else if (type == "CO")
	{
		p = make_shared<Particle_CO>();
	}
	else if (type == "CO2")
	{
		p = make_shared<Particle_CO2>();
	}
	else
	{
		cout << "Invalid particle type specified! Please check configuration file.\n";
		exit(1);
	}
	return p;
}
