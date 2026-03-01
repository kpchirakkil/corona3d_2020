/*
 * Background_Species.cpp
 *
 *  Created on: Jun 29, 2020
 *      Author: rodney
 */

#include "Background_Species.hpp"
#include <cctype>
#include <iomanip>
#include <utility>

namespace {
	bool is_non_finite_token(const string &token)
	{
		string lower = token;
		for (auto &c : lower)
		{
			c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
		}

		return (lower == "nan" || lower == "+nan" || lower == "-nan" ||
		        lower == "inf" || lower == "+inf" || lower == "-inf" ||
		        lower == "infinity" || lower == "+infinity" || lower == "-infinity");
	}

	bool parse_numeric_token(const string &token, double &value)
	{
		size_t pos = 0;
		try
		{
			value = stod(token, &pos);
		}
		catch (...)
		{
			return false;
		}

		return (pos > 0);
	}

	bool is_integer_token(const string &token)
	{
		if (token.empty())
		{
			return false;
		}

		size_t start = 0;
		if (token[0] == '+' || token[0] == '-')
		{
			start = 1;
		}
		if (start >= token.size())
		{
			return false;
		}

		for (size_t i = start; i < token.size(); i++)
		{
			if (!std::isdigit(static_cast<unsigned char>(token[i])))
			{
				return false;
			}
		}
		return true;
	}

	bool build_scattering_cdf(const vector<double> &angle_deg, const vector<double> &sigma,
	                          vector<double> &cdf, vector<double> &theta_rad)
	{
		if (angle_deg.empty() || angle_deg.size() != sigma.size())
		{
			return false;
		}

		size_t num_angles = angle_deg.size();
		cdf.resize(num_angles);
		theta_rad.resize(num_angles);
		vector<double> weighted_sigma(num_angles, 0.0);

		double sig_total = 0.0;
		for (size_t i = 0; i < num_angles; i++)
		{
			theta_rad[i] = angle_deg[i] * (constants::pi / 180.0);
			weighted_sigma[i] = sigma[i] * sin(theta_rad[i]);
			sig_total += weighted_sigma[i];
		}

		if (sig_total <= 0.0)
		{
			return false;
		}

		double running = 0.0;
		for (size_t i = 0; i < num_angles; i++)
		{
			running += weighted_sigma[i] / sig_total;
			cdf[i] = running;
		}
		cdf.back() = 1.0;
		return true;
	}

	bool load_inelastic_channel_table(const string &filename, vector<InelasticChannel> &channels)
	{
		ifstream infile;
		infile.open(filename);
		if (!infile.good())
		{
			return false;
		}

		channels.clear();
		size_t line_number = 0;
		string line;
		while (getline(infile, line))
		{
			line_number++;
			if (line.empty() || line[0] == '#' ||
			    std::all_of(line.begin(), line.end(), [](unsigned char c){ return std::isspace(c); }))
			{
				continue;
			}

			replace(line.begin(), line.end(), ',', ' ');
			stringstream str(line);

			string t_ji;
			if (!(str >> t_ji))
			{
				continue;
			}

			// Allow a header row like "ji,jf,sigma_cm2,delta_E_eV".
			double ji_val = 0.0;
			if (!parse_numeric_token(t_ji, ji_val) || is_non_finite_token(t_ji))
			{
				continue;
			}

			string t_jf, t_sigma, t_delta;
			if (!(str >> t_jf >> t_sigma >> t_delta))
			{
				cout << "ERROR: Could not parse inelastic channel table \"" << filename
				     << "\" at line " << line_number << ".\n";
				exit(1);
			}

			double jf_val = 0.0;
			double sigma_val = 0.0;
			double delta_val = 0.0;
			if (!parse_numeric_token(t_jf, jf_val) || is_non_finite_token(t_jf) ||
			    !parse_numeric_token(t_sigma, sigma_val) || is_non_finite_token(t_sigma) ||
			    !parse_numeric_token(t_delta, delta_val) || is_non_finite_token(t_delta))
			{
				cout << "ERROR: Non-numeric value in inelastic channel table \"" << filename
				     << "\" at line " << line_number << ".\n";
				exit(1);
			}

			InelasticChannel channel;
			channel.ji = static_cast<int>(ji_val);
			channel.jf = static_cast<int>(jf_val);
			channel.sigma_cm2 = sigma_val;
			channel.delta_E_eV = delta_val;

			// Keep only state-changing channels in this table.
			if (channel.jf != channel.ji && channel.sigma_cm2 > 0.0)
			{
				channels.push_back(channel);
			}
		}

		infile.close();
		return !channels.empty();
	}

	bool load_inelastic_channel_angle_cdfs(const string &filename, unordered_map<int, InelasticChannelAngleCDF> &angle_cdfs)
	{
		ifstream infile;
		infile.open(filename);
		if (!infile.good())
		{
			return false;
		}

		angle_cdfs.clear();
		int current_ji = -1;
		int current_jf = -1;
		vector<double> angle_deg;
		vector<double> dcs;

		auto flush_channel = [&]() {
			if (current_ji == 0 && current_jf >= 0 && !angle_deg.empty())
			{
				InelasticChannelAngleCDF channel_cdf;
				if (build_scattering_cdf(angle_deg, dcs, channel_cdf.cdf, channel_cdf.theta_rad))
				{
					angle_cdfs[current_jf] = std::move(channel_cdf);
				}
			}
			angle_deg.clear();
			dcs.clear();
		};

		size_t line_number = 0;
		string line;
		while (getline(infile, line))
		{
			line_number++;
			if (line.empty() || line[0] == '#' ||
			    std::all_of(line.begin(), line.end(), [](unsigned char c){ return std::isspace(c); }))
			{
				continue;
			}

			stringstream str(line);
			vector<string> tokens;
			string token;
			while (str >> token)
			{
				tokens.push_back(token);
			}
			if (tokens.size() < 3)
			{
				continue;
			}

			// Block header format: E ji jf
			if (is_integer_token(tokens[1]) && is_integer_token(tokens[2]))
			{
				flush_channel();
				current_ji = stoi(tokens[1]);
				current_jf = stoi(tokens[2]);
				continue;
			}

			if (current_ji < 0 || current_jf < 0)
			{
				continue;
			}

			double theta = 0.0;
			double sigma_val = 0.0;
			if (!parse_numeric_token(tokens[0], theta) || is_non_finite_token(tokens[0]) ||
			    !parse_numeric_token(tokens[1], sigma_val) || is_non_finite_token(tokens[1]))
			{
				cout << "ERROR: Could not parse DCS line in \"" << filename
				     << "\" at line " << line_number << ".\n";
				exit(1);
			}

			angle_deg.push_back(theta);
			dcs.push_back(sigma_val);
		}

		flush_channel();
		infile.close();
		return !angle_cdfs.empty();
	}

	string build_channel_table_path(const string &inelastic_dcs_path, int energy_index)
	{
		size_t slash = inelastic_dcs_path.find_last_of("/\\");
		string dir = (slash == string::npos) ? "" : inelastic_dcs_path.substr(0, slash + 1);
		ostringstream name;
		name << dir << "inelastic_channels_iEng" << setfill('0') << setw(2) << (energy_index + 1) << ".csv";
		return name.str();
	}

	string build_raw_inelastic_dcs_path(const string &inelastic_dcs_path, int energy_index)
	{
		size_t slash = inelastic_dcs_path.find_last_of("/\\");
		string inelastic_dir = (slash == string::npos) ? "" : inelastic_dcs_path.substr(0, slash + 1);  // .../inelastic/

		string inelastic_dir_trimmed = inelastic_dir;
		if (!inelastic_dir_trimmed.empty() &&
		    (inelastic_dir_trimmed.back() == '/' || inelastic_dir_trimmed.back() == '\\'))
		{
			inelastic_dir_trimmed.pop_back();
		}
		size_t parent_slash = inelastic_dir_trimmed.find_last_of("/\\");
		string base_dir = (parent_slash == string::npos) ? "" : inelastic_dir_trimmed.substr(0, parent_slash + 1); // .../O-CO2_full/

		ostringstream name;
		name << base_dir << "DCS_inelastic_3pes/DCS-allj_Eng_" << setfill('0') << setw(3) << (energy_index + 1) << ".dat";
		return name.str();
	}

	vector<double> get_last_finite_density_altitudes(const string &filename, int num_species)
	{
		ifstream infile;
		infile.open(filename);
		if (!infile.good())
		{
			cout << "\"" << filename << "\" not found!\n";
			exit(1);
		}

		vector<double> last_finite_alt(num_species, -1.0);
		size_t line_number = 0;
		string line;
		while (getline(infile, line))
		{
			line_number++;
			if (line.empty() || line[0] == '#' ||
			    std::all_of(line.begin(), line.end(), [](unsigned char c){ return std::isspace(c); }))
			{
				continue;
			}

			replace(line.begin(), line.end(), ',', ' ');
			stringstream str(line);
			string first_token;
			if (!(str >> first_token))
			{
				continue;
			}

			double alt = 0.0;
			if (is_non_finite_token(first_token))
			{
				cout << "ERROR: Non-finite altitude in \"" << filename << "\" at line " << line_number << ".\n";
				exit(1);
			}
			else if (!parse_numeric_token(first_token, alt))
			{
				continue;  // header row
			}

			for (int i=0; i<num_species; i++)
			{
				string token;
				if (!(str >> token))
				{
					cout << "ERROR: Could not parse column " << (i+2) << " in \"" << filename
					     << "\" at line " << line_number << ".\n";
					exit(1);
				}

				if (!is_non_finite_token(token))
				{
					double dummy = 0.0;
					if (!parse_numeric_token(token, dummy))
					{
						cout << "ERROR: Could not parse column " << (i+2) << " in \"" << filename
						     << "\" at line " << line_number << ".\n";
						exit(1);
					}
					last_finite_alt[i] = alt;
				}
			}
		}
		infile.close();

		for (int i=0; i<num_species; i++)
		{
			if (last_finite_alt[i] < 0.0)
			{
				cout << "ERROR: Density column " << (i+2) << " in \"" << filename
				     << "\" contains no finite values.\n";
				exit(1);
			}
		}

		return last_finite_alt;
	}
}

Background_Species::Background_Species() {
	use_temp_profile = false;
	use_dens_profile = false;
	profile_bottom_alt = 0.0;
	profile_top_alt = 0.0;
	num_collisions = 0;
	num_inelastic_collisions = 0;
	num_superelastic_collisions = 0;
	num_species = 0;
	ref_temp = 0.0;
	ref_height = 0.0;
	ref_g = 0.0;
	collision_target = -1;
	collision_theta = 0.0;
	my_dist = NULL;
	last_outcome = {false, false, -1, 0.0, 0.0, -1, -1};
}

Background_Species::Background_Species(int num_parts, string config_files[], Planet p, double ref_T, double ref_h, string temp_profile_filename, string dens_profile_filename, double profile_bottom, double profile_top)
{
	num_collisions = 0;
	num_inelastic_collisions = 0;
	num_superelastic_collisions = 0;
	last_outcome = {false, false, -1, 0.0, 0.0, -1, -1};
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
	dens_effective_top_alt.resize(num_species, profile_top_alt);
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
	missing_deltaE_use_avg.resize(num_species, true);
	inelastic_CDFs.resize(num_species);
	inelastic_channels.resize(num_species);
	inelastic_channel_angle_cdfs.resize(num_species);
	inelastic_rot_const_eV.resize(num_species, 0.0);
	inelastic_rot_pop_model.resize(num_species, InelasticRotPopModel::GroundStateJi0);
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
			else if (parameters[j] == "rot_const_cm1")
			{
				// Optional rotational constant for target internal-state population model.
				// Convert from cm^-1 to eV.
				inelastic_rot_const_eV[i] = stod(values[j]) * 1.2398419843320026e-4;
			}
			else if (parameters[j] == "rot_population_model")
			{
				if (values[j] == "ji0" || values[j] == "ground_state" || values[j] == "ground_state_ji0")
				{
					inelastic_rot_pop_model[i] = InelasticRotPopModel::GroundStateJi0;
				}
				else if (values[j] == "thermal" || values[j] == "thermal_boltzmann")
				{
					inelastic_rot_pop_model[i] = InelasticRotPopModel::ThermalBoltzmann;
				}
				else
				{
					cout << "ERROR: Unknown rot_population_model \"" << values[j]
					     << "\" for species config " << config_files[i]
					     << ". Use ji0 or thermal." << endl;
					exit(1);
				}
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
			else if (parameters[j] == "missing_deltaE_use_avg")
			{
				missing_deltaE_use_avg[i] = (values[j] == "true");
			}
		}

		// Default CO2 rotational constant if inelastic is enabled and no value is specified.
		if (enable_inelastic[i] && inelastic_rot_const_eV[i] <= 0.0 && bg_parts[i] && bg_parts[i]->get_name() == "CO2")
		{
			inelastic_rot_const_eV[i] = 0.39021 * 1.2398419843320026e-4;  // CO2 rotational constant (eV)
		}

		// Load inelastic DCS files if inelastic is enabled
		if (enable_inelastic[i])
		{
			inelastic_CDFs[i].resize(num_energies);
			inelastic_channels[i].resize(num_energies);
			inelastic_channel_angle_cdfs[i].resize(num_energies);
			int loaded_channel_tables = 0;
			int loaded_channel_angle_tables = 0;
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

						// Optional state-resolved inelastic channel table.
						string channel_path = build_channel_table_path(values[j], eidx);
						vector<InelasticChannel> channels;
						if (load_inelastic_channel_table(channel_path, channels))
						{
							inelastic_channels[i][eidx] = std::move(channels);
							loaded_channel_tables++;
						}

						// Optional channel-resolved inelastic angular DCS tables.
						// Currently these are available for ji=0 channels in DCS-allj files.
						string raw_dcs_path = build_raw_inelastic_dcs_path(values[j], eidx);
						unordered_map<int, InelasticChannelAngleCDF> channel_angle_cdfs;
						if (load_inelastic_channel_angle_cdfs(raw_dcs_path, channel_angle_cdfs))
						{
							inelastic_channel_angle_cdfs[i][eidx] = std::move(channel_angle_cdfs);
							loaded_channel_angle_tables++;
						}
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
				if (loaded_channel_tables > 0)
				{
					cout << "State-resolved inelastic channel tables loaded for species "
					     << bg_parts[i]->get_name() << " (" << loaded_channel_tables
					     << "/" << num_energies << " energies)." << endl;
				}
				if (loaded_channel_angle_tables > 0)
				{
					cout << "State-resolved inelastic angle tables (ji=0 channels) loaded for species "
					     << bg_parts[i]->get_name() << " (" << loaded_channel_angle_tables
					     << "/" << num_energies << " energies)." << endl;
				}
				if (missing_deltaE_use_avg[i])
				{
					cout << "Missing ji>0 channel delta_E fallback enabled for species "
					     << bg_parts[i]->get_name() << " (uses avg_energy_loss_file when sampled delta_E is 0)." << endl;
				}
				else
				{
					cout << "Missing ji>0 channel delta_E fallback disabled for species "
					     << bg_parts[i]->get_name() << " (strict table-only delta_E)." << endl;
				}
				if (missing_deltaE_use_avg[i] && !avg_eloss_interp[i])
				{
					cout << "WARNING: missing_deltaE_use_avg=true but avg_energy_loss_file is not loaded for species "
					     << bg_parts[i]->get_name() << ". Missing ji>0 channel delta_E values remain 0." << endl;
				}
				cout << "Rotational population model for species " << bg_parts[i]->get_name() << ": "
				     << (inelastic_rot_pop_model[i] == InelasticRotPopModel::GroundStateJi0 ?
				         "ji=0 only" : "thermal Boltzmann") << endl;
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
		vector<double> last_finite_alt = get_last_finite_density_altitudes(dens_profile_filename, num_species);
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
		for (int i=0; i<num_species; i++)
		{
			dens_effective_top_alt[i] = min(profile_top_alt, last_finite_alt[i]);
			if (dens_effective_top_alt[i] < profile_bottom_alt)
			{
				cout << "ERROR: Effective top altitude for species " << bg_parts[i]->get_name()
				     << " is below profile bottom altitude.\n";
				exit(1);
			}

			// generate interpolator for each density profile
			dens_interp[i] = make_shared<Interpolator>(dens_alt_bins, bg_densities[i]);

			// calc top and bottom scale height to be used for extrapolating densities
			bg_scaleheights[i][0] = constants::k_b*Tn_interp->loglinterp(profile_bottom_alt)/(bg_parts[i]->get_mass()*bottom_local_g);
			double top_local_g = (constants::G * my_planet.get_mass()) /
			                    (pow(my_planet.get_radius()+dens_effective_top_alt[i], 2.0));
			bg_scaleheights[i][1] = constants::k_b*Tn_interp->loglinterp(dens_effective_top_alt[i]) /
			                    (bg_parts[i]->get_mass()*top_local_g);

			if (dens_effective_top_alt[i] < profile_top_alt)
			{
				cout << "WARNING: Using effective density top for " << bg_parts[i]->get_name()
				     << " at " << dens_effective_top_alt[i]/1e5 << " km (profile top = "
				     << profile_top_alt/1e5 << " km)." << endl;
			}
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
			last_outcome.ji = -1;
			last_outcome.jf = -1;

			if (enable_inelastic[collision_target] && elastic_frac_interp[collision_target])
			{
				double f_el = elastic_frac_interp[collision_target]->linterp(energy[collision_target]);
				f_el = max(0.0, min(1.0, f_el));

					if (common::get_rand() > f_el)
					{
						// INELASTIC
						last_outcome.is_inelastic = true;
						bool sampled = sample_inelastic_transition(
							collision_target,
							energy[collision_target],
							alt,
							last_outcome.theta,
							last_outcome.delta_E_eV,
							last_outcome.ji,
							last_outcome.jf);

						// Legacy fallback (effectively dead code when rot_B > 0):
						// delta_E is now computed from quantum numbers B*[jf(jf+1)-ji(ji+1)]
						// in sample_inelastic_transition(), so delta_E_eV == 0 only when
						// ji == jf (filtered out) or rot_B == 0 (no rotational constant).
						if (sampled &&
						    collision_target >= 0 && collision_target < (int)missing_deltaE_use_avg.size() &&
						    missing_deltaE_use_avg[collision_target] &&
						    last_outcome.ji > 0 &&
						    last_outcome.delta_E_eV == 0.0 &&
						    avg_eloss_interp[collision_target])
						{
							last_outcome.delta_E_eV = avg_eloss_interp[collision_target]->linterp(energy[collision_target]);
						}

						// Fallback to legacy averaged model if channel-resolved data is unavailable.
						// avg_energy_loss values are also legacy (will use rotational formula in future).
						if (!sampled)
						{
							if (avg_eloss_interp[collision_target])
							{
								last_outcome.delta_E_eV = avg_eloss_interp[collision_target]->linterp(energy[collision_target]);
							}
							else
							{
								last_outcome.delta_E_eV = 0.0;
							}
							last_outcome.theta = find_new_theta_inelastic(collision_target, energy[collision_target]);
						}

					if (last_outcome.delta_E_eV < 0.0)
					{
						num_superelastic_collisions++;
					}
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
			last_outcome = {false, false, -1, 0.0, 0.0, -1, -1};
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

double Background_Species::get_local_neutral_temp(double alt)
{
	if (!use_temp_profile || !Tn_interp || Tn.empty())
	{
		return ref_temp;
	}

	if (alt < profile_bottom_alt)
	{
		return Tn.front();
	}
	if (alt > profile_top_alt)
	{
		return Tn.back();
	}
	return Tn_interp->loglinterp(alt);
}

bool Background_Species::sample_inelastic_transition(int part_index, double energy, double alt, double &theta, double &delta_E_eV, int &ji, int &jf)
{
	theta = 0.0;
	delta_E_eV = 0.0;
	ji = -1;
	jf = -1;

	if (part_index < 0 || part_index >= (int)diff_sigma_energies.size())
	{
		return false;
	}

	int num_energies = diff_sigma_energies[part_index].size();
	if (num_energies == 0)
	{
		return false;
	}

	int energy_index = 0;
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

	if (part_index >= (int)inelastic_channels.size() ||
	    energy_index >= (int)inelastic_channels[part_index].size())
	{
		return false;
	}

	const vector<InelasticChannel> &channels = inelastic_channels[part_index][energy_index];
	if (channels.empty())
	{
		return false;
	}

	int max_ji = 0;
	for (const auto &ch : channels)
	{
		if (ch.ji > max_ji)
		{
			max_ji = ch.ji;
		}
	}
	if (max_ji < 0)
	{
		return false;
	}

	vector<double> state_pop(max_ji + 1, 0.0);
	const double rot_B_eV = (part_index < (int)inelastic_rot_const_eV.size()) ? inelastic_rot_const_eV[part_index] : 0.0;
	InelasticRotPopModel pop_model = InelasticRotPopModel::GroundStateJi0;
	if (part_index >= 0 && part_index < (int)inelastic_rot_pop_model.size())
	{
		pop_model = inelastic_rot_pop_model[part_index];
	}

	if (pop_model == InelasticRotPopModel::GroundStateJi0)
	{
		state_pop[0] = 1.0;
	}
	else
	{
		const double local_T = max(1.0, get_local_neutral_temp(alt));
		const double k_B_eV = 8.617333262145e-5;

		if (rot_B_eV > 0.0)
		{
			double pop_sum = 0.0;
			for (int j=0; j<=max_ji; j++)
			{
				double E_j = rot_B_eV * static_cast<double>(j) * (static_cast<double>(j) + 1.0);
				double weight = (2.0 * static_cast<double>(j) + 1.0) * exp(-E_j / (k_B_eV * local_T));
				state_pop[j] = weight;
				pop_sum += weight;
			}

			if (pop_sum > 0.0)
			{
				for (int j=0; j<=max_ji; j++)
				{
					state_pop[j] /= pop_sum;
				}
			}
			else
			{
				double uniform = 1.0 / static_cast<double>(max_ji + 1);
				for (int j=0; j<=max_ji; j++)
				{
					state_pop[j] = uniform;
				}
			}
		}
		else
		{
			double uniform = 1.0 / static_cast<double>(max_ji + 1);
			for (int j=0; j<=max_ji; j++)
			{
				state_pop[j] = uniform;
			}
		}
	}

	vector<int> active_channel_indices;
	vector<double> channel_cdf;
	active_channel_indices.reserve(channels.size());
	channel_cdf.reserve(channels.size());

	double total_weight = 0.0;
	for (int idx=0; idx<(int)channels.size(); idx++)
	{
		const InelasticChannel &ch = channels[idx];
		if (ch.jf == ch.ji || ch.sigma_cm2 <= 0.0)
		{
			continue;
		}

		if (ch.ji < 0 || ch.ji > max_ji)
		{
			continue;
		}

		double weight = ch.sigma_cm2 * state_pop[ch.ji];
		if (weight <= 0.0)
		{
			continue;
		}

		total_weight += weight;
		active_channel_indices.push_back(idx);
		channel_cdf.push_back(total_weight);
	}

	if (total_weight <= 0.0 || channel_cdf.empty())
	{
		return false;
	}

	double u = common::get_rand() * total_weight;
	int picked = 0;
	while (picked + 1 < (int)channel_cdf.size() && u > channel_cdf[picked])
	{
		picked++;
	}

	const InelasticChannel &chosen = channels[active_channel_indices[picked]];
	ji = chosen.ji;
	jf = chosen.jf;
	// Compute delta_E from rotational quantum numbers: B * [jf(jf+1) - ji(ji+1)]
	// This is the correct internal energy change; the CSV delta_E_eV values are legacy.
	if (rot_B_eV > 0.0)
	{
		delta_E_eV = rot_B_eV * (static_cast<double>(jf) * (static_cast<double>(jf) + 1.0)
		                        - static_cast<double>(ji) * (static_cast<double>(ji) + 1.0));
	}
	else
	{
		delta_E_eV = chosen.delta_E_eV;  // fallback if no rotational constant
	}

	// Hybrid angle sampling:
	// - For ji=0: use channel-resolved angle CDF keyed by jf
	// - For ji>0: use ji=0 DCS as proxy, keyed by jf first, then |jf-ji| (same delta_j)
	// - Fallback: species/energy aggregate inelastic DCS
	bool sampled_channel_angle = false;
	if (part_index >= 0 && part_index < (int)inelastic_channel_angle_cdfs.size() &&
	    energy_index >= 0 && energy_index < (int)inelastic_channel_angle_cdfs[part_index].size())
	{
		auto &energy_angle_tables = inelastic_channel_angle_cdfs[part_index][energy_index];

		// Try direct lookup by jf (works for ji=0; for ji>0 may match if same jf exists)
		auto it = energy_angle_tables.find(jf);

		// For ji>0, if direct jf lookup fails, try |jf - ji| as proxy key
		// (ji=0 channel with same magnitude of angular momentum change)
		if (it == energy_angle_tables.end() && ji > 0)
		{
			int delta_j = std::abs(jf - ji);
			it = energy_angle_tables.find(delta_j);
		}

		if (it != energy_angle_tables.end() && !it->second.cdf.empty() && !it->second.theta_rad.empty())
		{
			double u_theta = common::get_rand();
			u_theta = std::max(0.0, std::min(u_theta, 1.0));
			int k = 0;
			while ((k + 1 < (int)it->second.cdf.size()) && it->second.cdf[k] < u_theta)
			{
				k++;
			}
			theta = it->second.theta_rad[k];
			sampled_channel_angle = true;
		}
	}
	if (!sampled_channel_angle)
	{
		theta = find_new_theta_inelastic(part_index, energy);
	}

	return true;
}

// get density from imported density profile
double Background_Species::get_density(double alt, int index)
{
	double current_dens = 0.0;
	double top_alt = profile_top_alt;
	if (index >= 0 && index < (int)dens_effective_top_alt.size())
	{
		top_alt = dens_effective_top_alt[index];
	}

	// if outside of profile boundaries, need to extrapolate using a scale height
	if (alt < profile_bottom_alt)
	{
		current_dens = calc_new_density(bg_densities[index][0], bg_scaleheights[index][0], profile_bottom_alt - alt);
	}
	else if (alt > top_alt)
	{
		double dens_at_top = dens_interp[index]->loglinterp(top_alt);
		current_dens = calc_new_density(dens_at_top, bg_scaleheights[index][1], top_alt - alt);
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

int Background_Species::get_num_superelastic_collisions()
{
	return num_superelastic_collisions;
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
