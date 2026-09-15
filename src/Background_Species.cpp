/*
 * Background_Species.cpp
 *
 *  Created on: Jun 29, 2020
 *      Author: rodney
 */

#include "Background_Species.hpp"
#include <cctype>
#include <limits>
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

    double angular_integral(const vector<double> &angle, const vector<double> &sigma) {
        double sum=0, prev_theta=0, prev_weight=0;
        for (size_t k=0;k<angle.size();++k) {
            const double t=angle[k]*constants::pi/180;
            const double w=(angle[k]==180)?0.0:sigma[k]*sin(t);
            sum+=0.5*(t-prev_theta)*(w+prev_weight);
            prev_theta=t; prev_weight=w;
        }
        sum+=0.5*(constants::pi-prev_theta)*prev_weight;
        return constants::twopi*sum;
    }

    bool build_scattering_cdf(const vector<double> &angle_deg, const vector<double> &sigma,
                              vector<double> &cdf, vector<double> &theta_rad) {
        if (angle_deg.empty() || angle_deg.size()!=sigma.size()) return false;
        cdf.clear(); theta_rad.clear();
        double prev_theta=0,prev_weight=0,sum=0;
        cdf.push_back(0); theta_rad.push_back(0);
        for (size_t k=0;k<angle_deg.size();++k) {
            if (!std::isfinite(sigma[k]) || sigma[k]<0 || angle_deg[k]<0 || angle_deg[k]>180 ||
                (k && angle_deg[k]<=angle_deg[k-1]))
                throw std::runtime_error("Invalid angular cross-section grid");
            double t=angle_deg[k]*constants::pi/180;
            double w=(angle_deg[k]==180)?0.0:sigma[k]*sin(t);
            sum+=0.5*(t-prev_theta)*(w+prev_weight);
            if (t>prev_theta) { cdf.push_back(sum); theta_rad.push_back(t); }
            prev_theta=t; prev_weight=w;
        }
        if (prev_theta<constants::pi) {
            sum+=0.5*(constants::pi-prev_theta)*prev_weight;
            cdf.push_back(sum); theta_rad.push_back(constants::pi);
        }
        if (!(sum>0) || !std::isfinite(sum)) return false;
        for (double &x:cdf) x/=sum;
        cdf.back()=1;
        return true;
    }

    double draw_angle(const vector<double> &cdf, const vector<double> &angles) {
        if (cdf.empty()) throw std::runtime_error("Empty scattering CDF");
        const double u=common::get_rand();
        size_t hi=std::upper_bound(cdf.begin(),cdf.end(),u)-cdf.begin();
        if (hi>=cdf.size()) return angles.back();
        if (hi==0) return angles.front();
        const double t=(u-cdf[hi-1])/(cdf[hi]-cdf[hi-1]);
        return angles[hi-1]+t*(angles[hi]-angles[hi-1]);
    }

    int nearest_energy(const vector<double> &grid, double e) {
        auto hi=std::lower_bound(grid.begin(),grid.end(),e);
        if (hi==grid.begin()) return 0;
        if (hi==grid.end()) return grid.size()-1;
        return (e-*(hi-1)<=*hi-e)?int(hi-grid.begin()-1):int(hi-grid.begin());
    }

    bool load_inelastic_channel_angle_cdfs(const string &filename, unordered_map<int, InelasticChannelAngleCDF> &angle_cdfs)
    {
        ifstream input(filename);
        if (!input) throw std::runtime_error("Cannot open configured channel-angle file: " + filename);
        angle_cdfs.clear();
        int ji=-1, jf=-1;
        vector<double> angles, sigma;
        unordered_map<int,bool> seen;
        size_t line_number=0;
        auto invalid=[&](const string &reason) {
            return std::runtime_error("Invalid channel-angle file " + filename + " at line " +
                                      to_string(line_number) + ": " + reason);
        };
        auto flush=[&]() {
            if (ji<0) return;
            if (angles.size()<2) throw invalid("block needs at least two angular rows");
            InelasticChannelAngleCDF cdf;
            // Zero-weight closed channels are allowed, but supply no angle CDF.
            if (build_scattering_cdf(angles,sigma,cdf.cdf,cdf.theta_rad) && jf!=ji)
                angle_cdfs.emplace(jf,std::move(cdf));
            angles.clear(); sigma.clear();
        };
        string line;
        while (getline(input,line)) {
            ++line_number;
            const auto first=line.find_first_not_of(" \t\r");
            if (first==string::npos || line[first]=='#') continue;
            istringstream row(line);
            string a,b,c,extra;
            if (!(row>>a>>b>>c) || (row>>extra)) throw invalid("expected exactly three columns");
            double x,y,z;
            auto numeric=[](const string &token,double &value) {
                size_t used=0;
                try { value=stod(token,&used); } catch (...) { return false; }
                return used==token.size() && std::isfinite(value);
            };
            if (!numeric(a,x) || !numeric(b,y) || !numeric(c,z)) throw invalid("nonfinite or nonnumeric value");
            // Block headers have integer ji/jf; angular rows use decimal or
            // scientific notation for DCS and Te, as in the source data.
            if (is_integer_token(b) && is_integer_token(c)) {
                flush();
                if (x<0 || y!=0 || z<0 || z>100000) throw invalid("only nonnegative-energy ji=0 blocks are supported");
                ji=0; jf=static_cast<int>(z);
                if (!seen.emplace(jf,true).second) throw invalid("duplicate final state");
            } else {
                if (ji<0) throw invalid("angular row before block header");
                if (x<0 || x>180 || y<0 || (!angles.empty() && x<=angles.back()))
                    throw invalid("invalid angle grid or negative DCS");
                angles.push_back(x); sigma.push_back(y);
                // Te is validated but not used as internal energy transfer.
                // collision_channels() computes B*[jf(jf+1)-ji(ji+1)].
            }
        }
        flush();
        return !angle_cdfs.empty();
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
	inelastic_CDFs.resize(num_species);
	inelastic_channel_angle_cdfs.resize(num_species);
	inelastic_rot_const_eV.resize(num_species, 0.0);
	inelastic_rot_pop_model.resize(num_species, InelasticRotPopModel::GroundStateJi0);
    rotational_tables.resize(num_species);
    collision_sigma_bounds.resize(num_species);
    elastic_dcs_integrals.resize(num_species);
    thermal_angle_proxy.resize(num_species,false);
    inelastic_angle_min_eV.resize(num_species,0.0);
	for (int i=0; i<num_species; i++)
	{
		string rotational_file;
        string inelastic_model="state_resolved";
        int population_max_j=-1;
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
            else if (parameters[j] == "inelastic_model") {
                inelastic_model=values[j];
                if (inelastic_model!="state_resolved" && inelastic_model!="legacy_average")
                    throw std::runtime_error("inelastic_model must be state_resolved or legacy_average: " + config_files[i]);
            }
            else if (parameters[j] == "rotational_cross_sections_file") rotational_file=values[j];
            else if (parameters[j] == "rot_population_max_j") {
                size_t used=0;
                population_max_j=stoi(values[j],&used);
                if (used!=values[j].size() || population_max_j<0)
                    throw std::runtime_error("rot_population_max_j must be a non-negative integer: " + config_files[i]);
            }
            else if (parameters[j] == "inelastic_angle_min_eV") inelastic_angle_min_eV[i]=stod(values[j]);
            else if (parameters[j] == "thermal_angular_model") {
                if (values[j]!="ji0_proxy") throw std::runtime_error("thermal_angular_model must be ji0_proxy");
                thermal_angle_proxy[i]=true;
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

		// Default CO2 rotational constant if inelastic is enabled and no value is specified.
		if (enable_inelastic[i] && inelastic_rot_const_eV[i] <= 0.0 && bg_parts[i] && bg_parts[i]->get_name() == "CO2")
		{
			inelastic_rot_const_eV[i] = 0.39021 * 1.2398419843320026e-4;  // CO2 rotational constant (eV)
		}

        if (enable_inelastic[i]) {
            if (inelastic_model=="state_resolved" && rotational_file.empty())
                throw std::runtime_error("state_resolved requires rotational_cross_sections_file: " + config_files[i]);
            if (inelastic_model=="legacy_average") {
                if (!rotational_file.empty())
                    throw std::runtime_error("legacy_average conflicts with rotational_cross_sections_file: " + config_files[i]);
                if (!sigma_total_interp[i] || !elastic_frac_interp[i] || !avg_eloss_interp[i])
                    throw std::runtime_error("legacy_average requires total_sigma_file_total, elastic_fraction_file and avg_energy_loss_file: " + config_files[i]);
            }
            cout << "Inelastic rate model for " << bg_parts[i]->get_name() << ": " << inelastic_model << endl;
        }
        if (enable_inelastic[i] && !rotational_file.empty()) {
            rotational_tables[i].load(rotational_file,inelastic_rot_const_eV[i]);
            cout << "Integral state cross sections: " << rotational_file << endl;
            if (bg_parts[i]->get_name()=="CO")
                cout << "WARNING: CO channels are limited to published jf<=30; omitted high-j rates are unknown." << endl;
        }
        if (enable_inelastic[i] && inelastic_rot_pop_model[i]==InelasticRotPopModel::ThermalBoltzmann) {
            if (rotational_tables[i].max_ji==0)
                throw std::runtime_error("Thermal rotations need excited-initial-state integral data: " + config_files[i]);
            if (!thermal_angle_proxy[i])
                throw std::runtime_error("Only ji=0 angular data are available. Set thermal_angular_model ji0_proxy to explicitly accept this approximation.");
        }
        if (population_max_j>=0) {
            if (!enable_inelastic[i] || inelastic_rot_pop_model[i]!=InelasticRotPopModel::ThermalBoltzmann)
                throw std::runtime_error("rot_population_max_j requires rot_population_model thermal: " + config_files[i]);
            if (population_max_j>rotational_tables[i].max_ji)
                throw std::runtime_error("rot_population_max_j exceeds tabulated initial states: " + config_files[i]);
            rotational_tables[i].population_max_j=population_max_j;
            // Fail at load if a kept level lacks data; any positive temperature exercises this check.
            rotational_tables[i].populations(1000.0,true,bg_parts[i]->get_name());
        }
        if (enable_inelastic[i] && rotational_tables[i].curves.empty())
            cout << "WARNING: legacy averaged inelastic model for " << bg_parts[i]->get_name()
                 << "; outside-grid values are endpoint approximations; mean-loss channels are closed below their mean energy." << endl;

		// Load inelastic DCS files if inelastic is enabled
		if (enable_inelastic[i])
		{
			inelastic_CDFs[i].resize(num_energies);
			inelastic_channel_angle_cdfs[i].resize(num_energies);
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
					}
				}
			}
            // Optional ji=0 channel angles use explicit paths and the documented
            // E ji jf / theta DCS Te block format, independent of species layout.
            const string angle_suffix="_inelastic_angles_file";
            for (int j=0; j<num_params; ++j) {
                const string &key=parameters[j];
                if (key.size()<=angle_suffix.size() ||
                    key.compare(key.size()-angle_suffix.size(),angle_suffix.size(),angle_suffix)!=0) continue;
                if (key.size()<=6+angle_suffix.size() || key.compare(0,6,"energy")!=0)
                    throw std::runtime_error("Invalid channel-angle key: " + key);
                const string index=key.substr(6,key.size()-6-angle_suffix.size());
                if (!is_integer_token(index))
                    throw std::runtime_error("Invalid channel-angle key: " + key);
                const int eidx=stoi(index)-1;
                if (eidx<0 || eidx>=num_energies)
                    throw std::runtime_error("Channel-angle energy index out of range: " + key);
                auto &angles=inelastic_channel_angle_cdfs[i][eidx];
                if (!angles.empty()) throw std::runtime_error("Duplicate channel-angle key: " + key);
                if (!load_inelastic_channel_angle_cdfs(values[j],angles))
                    throw std::runtime_error("Configured channel-angle file has no usable ji=0 inelastic angles: " + values[j]);
                ++loaded_channel_angle_tables;
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
				if (loaded_channel_angle_tables > 0)
				{
					cout << "State-resolved inelastic angle tables (ji=0 channels) loaded for species "
					     << bg_parts[i]->get_name() << " (" << loaded_channel_angle_tables
					     << "/" << num_energies << " energies)." << endl;
				}
                cout << "Inelastic angular model for " << bg_parts[i]->get_name() << ": "
                     << (loaded_channel_angle_tables>0?"configured ji=0 transition angles; aggregate fallback for uncovered channels":"aggregate angles for all transitions") << endl;
                if (inelastic_rot_pop_model[i]==InelasticRotPopModel::ThermalBoltzmann)
                    cout << "Thermal angles: reciprocal ji0_proxy at total energy; 0<->j channel DCS when available, aggregate otherwise (explicitly accepted)" << endl;
                cout << "Rotational population model: "
                     << (inelastic_rot_pop_model[i]==InelasticRotPopModel::GroundStateJi0?"ji=0":"thermal");
                if (rotational_tables[i].population_max_j>=0)
                    cout << " (Boltzmann truncated to j<=" << rotational_tables[i].population_max_j << " and renormalized)";
                cout << endl;
        }

		bg_scaleheights[i].push_back(constants::k_b*ref_temp/(bg_parts[i]->get_mass()*ref_g));
		bg_avg_v[i].push_back(sqrt(constants::k_b*ref_temp/bg_parts[i]->get_mass()));

		elastic_dcs_integrals[i].resize(num_energies);
		for (int j=0; j<num_energies; j++)
		{
			diff_sigma_CDFs[i][j].resize(2);
			vector<vector<double>> diff_sigma_PDF;
			diff_sigma_PDF.resize(2);
			diff_sigma_energies[i].push_back(stod(values[energies_index + j]));
			common::import_csv(values[energies_index + num_energies + j], diff_sigma_PDF[0], diff_sigma_PDF[1]);
            elastic_dcs_integrals[i][j]=angular_integral(diff_sigma_PDF[0],diff_sigma_PDF[1]);
			make_new_CDF(i, j, diff_sigma_PDF[0], diff_sigma_PDF[1]);
		}
        double elastic_max=bg_sigma_defaults[i];
        if (sigma_interp[i]) elastic_max=*std::max_element(bg_sigma_tables[i][1].begin(),bg_sigma_tables[i][1].end());
        if (!rotational_tables[i].curves.empty()) collision_sigma_bounds[i]=rotational_tables[i].upper_bound(elastic_max);
        else {
            double inelastic_max=0;
            if (enable_inelastic[i] && sigma_total_interp[i]) inelastic_max=sigma_total_interp[i]->max_value();
            collision_sigma_bounds[i]=(elastic_max+inelastic_max)*(1+1e-12);
        }
        if (!(collision_sigma_bounds[i]>0) || !std::isfinite(collision_sigma_bounds[i]))
            throw std::runtime_error("Invalid collision cross-section bound");
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
    const double m=p1->get_mass(), M=p2->get_mass();
    const Vector3d g(p1->get_vx()-p2->get_vx(),p1->get_vy()-p2->get_vy(),p1->get_vz()-p2->get_vz());
    return 0.5*(m*M/(m+M))*g.squaredNorm()/constants::ergev;
}

// calculates new density of background particle based on radial position and scale height
double Background_Species::calc_new_density(double ref_density, double scale_height, double r_moved)
{
	return ref_density*exp(r_moved/scale_height);
}

// State weights used by both the total collision rate and event selection.
vector<RotationalTransition> Background_Species::collision_channels(int i, double energy, double alt)
{
    double elastic=sigma_interp[i]?sigma_interp[i]->linterp(energy):bg_sigma_defaults[i];
    if (!std::isfinite(elastic) || elastic<0) throw std::runtime_error("Invalid elastic cross section");
    if (!enable_inelastic[i]) return {{-1,-1,elastic,0.0}};
    if (!rotational_tables[i].curves.empty())
        return rotational_tables[i].evaluate(energy,get_local_neutral_temp(alt),
            inelastic_rot_pop_model[i]==InelasticRotPopModel::ThermalBoltzmann,
            bg_parts[i]->get_name(),elastic);

    // Compatibility for older averaged configurations. Retain the independently
    // tabulated elastic rate; only the inelastic part comes from legacy scalars.
    vector<RotationalTransition> channels={{-1,-1,elastic,0.0}};
    if (sigma_total_interp[i] && elastic_frac_interp[i] && avg_eloss_interp[i]) {
        const double fraction=std::max(0.0,std::min(1.0,elastic_frac_interp[i]->linterp(energy)));
        const double inelastic=sigma_total_interp[i]->linterp(energy)*(1-fraction);
        const double de=avg_eloss_interp[i]->linterp(energy);
        if (inelastic>0 && de>0 && de<energy) channels.push_back({-1,-2,inelastic,de});
    }
    return channels;
}

// Draw the first physical event in dt using null-collision thinning.
// Proposal q(w) is proportional to f_MB(w)*(|v|+|w|), which bounds
// f_MB(w)*|v-w|. This samples the full relative-speed collision kernel,
// without a velocity cutoff or a finite-dt Bernoulli approximation.
bool Background_Species::check_collision(shared_ptr<Particle> p, double dt)
{
    last_outcome={false,false,-1,0.0,0.0,-1,-1};
    collision_target=-1;
    collision_delay=0;
    if (!std::isfinite(dt) || dt<0) throw std::invalid_argument("Invalid collision timestep");
    if (dt==0) return false;
    const double alt=p->get_radius()-my_planet.get_radius();
    const double v=p->get_total_v();
    const double temperature=get_local_neutral_temp(alt);
    if (!std::isfinite(temperature) || temperature<0) throw std::runtime_error("Invalid neutral temperature");
    vector<double> scales(num_species),means(num_species),rates(num_species);
    double total_rate=0;
    for (int i=0;i<num_species;++i) {
        double density=use_dens_profile?get_density(alt,i):
            calc_new_density(bg_densities[i][0],bg_scaleheights[i][0],ref_height-alt);
        if (!std::isfinite(density) || density<0) throw std::runtime_error("Invalid background density");
        scales[i]=sqrt(constants::k_b*temperature/bg_parts[i]->get_mass());
        means[i]=sqrt(8/constants::pi)*scales[i];
        rates[i]=density*collision_sigma_bounds[i]*(v+means[i]);
        total_rate+=rates[i];
    }
    if (!std::isfinite(total_rate)) throw std::runtime_error("Nonfinite collision rate");
    if (!(total_rate>0)) return false;
    while (true) {
        collision_delay+=-std::log1p(-common::get_rand())/total_rate;
        if (collision_delay>=dt) return false;
        double pick=common::get_rand()*total_rate;
        int i=0;
        while (i+1<num_species && pick>=rates[i]) { pick-=rates[i]; ++i; }
        if (common::get_rand()*(v+means[i])<v) {
            my_dist->init_vonly(bg_parts[i],scales[i]);
        } else {
            // Speed-biased Maxwell distribution: r^2/(2*s^2) ~ Gamma(2,1).
            const double speed=scales[i]*sqrt(-2*(std::log1p(-common::get_rand())+std::log1p(-common::get_rand())));
            const double z=2*common::get_rand()-1, phi=constants::twopi*common::get_rand();
            const double radial=speed*sqrt(std::max(0.0,1-z*z));
            bg_parts[i]->init_particle_vonly(radial*cos(phi),radial*sin(phi),speed*z);
        }
        const Vector3d relative(p->get_vx()-bg_parts[i]->get_vx(),p->get_vy()-bg_parts[i]->get_vy(),p->get_vz()-bg_parts[i]->get_vz());
        const double energy=calc_collision_e(p,bg_parts[i]);
        auto channels=collision_channels(i,energy,alt);
        double sigma=0;
        for (const auto &c:channels) sigma+=c.sigma;
        const double denominator=collision_sigma_bounds[i]*(v+bg_parts[i]->get_total_v());
        const double acceptance=denominator>0?relative.norm()*sigma/denominator:0;
        if (!std::isfinite(acceptance) || acceptance>1+1e-10 || acceptance<0)
            throw std::runtime_error("Invalid null-collision bound");
        if (common::get_rand()>=acceptance || !(sigma>0)) continue;
        pick=common::get_rand()*sigma;
        size_t k=0;
        while (k+1<channels.size() && pick>=channels[k].sigma) { pick-=channels[k].sigma; ++k; }
        const auto &chosen=channels[k];
        collision_target=i;
        const bool inelastic=chosen.ji!=chosen.jf;
        collision_theta=state_scattering_angle(i,energy,chosen.ji,chosen.jf);
        last_outcome={true,inelastic,i,collision_theta,chosen.delta_e,chosen.ji,chosen.jf};
        ++num_collisions;
        if (inelastic) ++num_inelastic_collisions;
        if (chosen.delta_e<0) ++num_superelastic_collisions;
        return true;
    }
}

void Background_Species::do_collisions(shared_ptr<Particle> p, double dt, double time)
{
    double elapsed=0;
    while (elapsed<dt && check_collision(p,dt-elapsed)) {
        elapsed+=collision_delay;
        const auto outcome=last_outcome;
        p->do_collision(get_collision_target(),outcome.theta,time+elapsed,my_planet.get_radius(),
                        outcome.is_inelastic,outcome.delta_E_eV,outcome.ji,outcome.jf);
    }
}

double Background_Species::find_new_theta(int i, double energy)
{
    const int k=nearest_energy(diff_sigma_energies[i],energy);
    // Preserve the resolved DCS contribution instead of scaling it to fill
    // missing forward scattering. The residual is represented as theta=0.
    // Its unresolved transport contribution is a documented angular-resolution
    // approximation; scalar elastic rates always use the integral reference.
    const double reference=sigma_interp[i]?sigma_interp[i]->linterp(diff_sigma_energies[i][k]):bg_sigma_defaults[i];
    const double resolved=elastic_dcs_integrals[i][k];
    if (reference>resolved && common::get_rand()*reference>=resolved) return 0.0;
    return draw_angle(diff_sigma_CDFs[i][k][0],diff_sigma_CDFs[i][k][1]);
}

double Background_Species::find_new_theta_inelastic(int i, double energy)
{
    energy=std::max(energy,inelastic_angle_min_eV[i]);
    const int k=nearest_energy(diff_sigma_energies[i],energy);
    if (k>=(int)inelastic_CDFs[i].size() || inelastic_CDFs[i][k].size()<2)
        throw std::runtime_error("Missing inelastic angular distribution");
    return draw_angle(inelastic_CDFs[i][k][0],inelastic_CDFs[i][k][1]);
}

double Background_Species::state_scattering_angle(int i, double energy, int ji, int jf)
{
    if (ji==jf) return find_new_theta(i,energy);
    // Legacy mean-loss events have no rotational states or total-energy grid.
    if (ji<0 || jf<0) return find_new_theta_inelastic(i,energy);

    // A reversible pair must share a normalized angular law at the same
    // total energy K + E(ji). Scalar cross sections already supply the
    // kinetic-energy and degeneracy factors required by detailed balance.
    energy+=rotational_tables[i].level(ji);
    // Reconstructing total energy in the reverse direction can differ by an
    // ulp. Canonicalize roundoff-sized midpoint ties before either lookup so
    // a pair cannot choose opposite sides of a discontinuous nearest bin.
    const auto &grid=diff_sigma_energies[i];
    const auto upper=std::lower_bound(grid.begin(),grid.end(),energy);
    if (upper!=grid.begin() && upper!=grid.end()) {
        const double midpoint=(*(upper-1)+*upper)/2;
        const double tolerance=8*std::numeric_limits<double>::epsilon()*std::max(energy,*upper);
        if (std::abs(energy-midpoint)<=tolerance) energy=midpoint;
    }
    const int low=std::min(ji,jf), high=std::max(ji,jf);
    if (low==0 && !inelastic_channel_angle_cdfs[i].empty()) {
        const int k=nearest_energy(diff_sigma_energies[i],energy);
        auto it=inelastic_channel_angle_cdfs[i][k].find(high);
        if (it!=inelastic_channel_angle_cdfs[i][k].end())
            return draw_angle(it->second.cdf,it->second.theta_rad);
    }
    // Missing ground-state DCS blocks and pairs with both states excited use
    // the same aggregate proxy in either direction, also at total energy.
    // This enforces reciprocity but cannot recover unavailable excited DCS.
    return find_new_theta_inelastic(i,energy);
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

void Background_Species::make_new_CDF(int i, int k, vector<double> &angle, vector<double> &sigma)
{
    if (!build_scattering_cdf(angle,sigma,diff_sigma_CDFs[i][k][0],diff_sigma_CDFs[i][k][1]))
        throw std::runtime_error("Zero or invalid elastic angular distribution");
}

void Background_Species::make_new_inelastic_CDF(int i, int k, vector<double> &angle, vector<double> &sigma)
{
    if (!build_scattering_cdf(angle,sigma,inelastic_CDFs[i][k][0],inelastic_CDFs[i][k][1]))
        throw std::runtime_error("Zero or invalid inelastic angular distribution");
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
