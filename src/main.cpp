/*
 * main.cpp
 *
 *  Created on: May 27, 2020
 *      Author: rodney
 */
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "constants.hpp"
#include "Atmosphere.hpp"
using namespace std;

//subroutine to set particle types from main
Particle* set_particle_type(string type)
{
	Particle* p;

	if (type == "H")
	{
		p = new Particle_H();
	}
	else if (type == "O")
	{
		p = new Particle_O();
	}
	else if (type == "N2")
	{
		p = new Particle_N2();
	}
	else if (type == "CO")
	{
		p = new Particle_CO();
	}
	else if (type == "CO2")
	{
		p = new Particle_CO2();
	}
	else
	{
		cout << "Invalid particle type specified! Please check configuration file.\n";
		exit(1);
	}
	return p;
}

int main(int argc, char* argv[])
{
	//initialize and read parameters from configuration file
	int num_testparts = 0;
	string part_type = "";
	string dist_type = "";
	string pos_infile = "";
	string vel_infile = "";
	int timesteps = 0;
	double dt = 0.0;
	double ref_height = 0.0;
	double bg_temp = 0.0;
	double planet_mass = 0.0;
	double planet_radius = 0.0;
	Planet my_planet;
	vector<Particle*> parts;
	Distribution* dist;
	int num_bgparts = 0;
	int bg_params_index = 0;

	ifstream infile;
	infile.open("/home/rodney/git/corona3d_2020/src/corona3d_2020.cfg");
	if (!infile.good())
	{
		cout << "Configuration file not found!\n";
		return 1;
	}
	string line, param, val;
	vector<string> parameters;
	vector<string> values;
	int num_params = 0;

	while (getline(infile, line))
	{
		if (line[0] == '#' || line.empty() || std::all_of(line.begin(), line.end(), ::isspace))
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

	for (int i=0; i<num_params; i++)
	{
		if (parameters[i] == "num_testparts")
		{
			num_testparts = stoi(values[i]);
		}
		else if (parameters[i] == "part_type")
		{
			part_type = values[i];
		}
		else if (parameters[i] == "dist_type")
		{
			dist_type = values[i];
		}
		else if (parameters[i] == "pos_infile")
		{
			pos_infile = values[i];
		}
		else if (parameters[i] == "vel_infle")
		{
			vel_infile = values[i];
		}
		else if (parameters[i] == "timesteps")
		{
			timesteps = stoi(values[i]);
		}
		else if (parameters[i] == "dt")
		{
			dt = stod(values[i]);
		}
		else if (parameters[i] == "ref_height")
		{
			ref_height = stod(values[i]);
		}
		else if (parameters[i] == "bg_temp")
		{
			bg_temp = stod(values[i]);
		}
		else if (parameters[i] == "planet_mass")
		{
			planet_mass = stod(values[i]);
		}
		else if (parameters[i] == "planet_radius")
		{
			planet_radius = stod(values[i]);
		}
		else if (parameters[i] == "num_bgparts")
		{
			num_bgparts = stoi(values[i]);
			bg_params_index = i+1;
		}
	}

	//initialize planet and test particles
	my_planet.init(planet_mass, planet_radius);
	parts.resize(num_testparts);
	for (int i=0; i<num_testparts; i++)
	{
		parts[i] = set_particle_type(part_type);
	}

	//instantiate the Distribution class to be used
	if (dist_type == "Hot_H")
	{
		//haven't made this distribution class yet
	}
	else if (dist_type == "Hot_O")
	{
		dist = new Distribution_Hot_O(my_planet, ref_height, bg_temp);
	}
	else if (dist_type == "MB")
	{
		dist = new Distribution_MB(my_planet, ref_height, bg_temp);
	}
	else if (dist_type == "Import")
	{
		dist = new Distribution_Import(my_planet, ref_height, bg_temp, pos_infile, vel_infile);
	}
	else
	{
		cout << "Invalid distribution type! Please check configuration file.\n";
		return 1;
	}

	//initialize Background_Species
	Particle* bg_parts[num_bgparts];
	double bg_dens[num_bgparts];
	double bg_sigs[num_bgparts];
	for (int i=0; i<num_bgparts; i++)
	{
		bg_parts[i] = set_particle_type(values[bg_params_index + i]);
		bg_dens[i] = stod(values[bg_params_index + num_bgparts + i]);
		bg_sigs[i] = stod(values[bg_params_index + num_bgparts + num_bgparts + i]);
	}
	Distribution_MB* bg_dist = new Distribution_MB(my_planet, ref_height, bg_temp);
	Background_Species bg_spec(num_bgparts, my_planet, bg_temp, ref_height, bg_dist, bg_parts, bg_dens, bg_sigs);

	// initialize atmosphere and run simulation
	Atmosphere my_atmosphere(num_testparts, my_planet, parts, dist, bg_spec, bg_temp, ref_height);
	my_atmosphere.output_velocity_distro(100.0, 150, "/home/rodney/Documents/coronaTest/vdist.out");
	my_atmosphere.run_simulation(dt, timesteps);
	my_atmosphere.output_velocity_distro(100.0, 150, "/home/rodney/Documents/coronaTest/vdist2.out");

/*
	// simulation parameters
	int num_testparts = 10000;
	int timesteps = 1000000;
	double dt = 0.05;
	double ref_height = 200e3;
	double bg_temp = 277.6;
	double planet_mass = 6.4185e23;
	double planet_radius = 3.397e6;
	Planet mars;
	mars.init(planet_mass, planet_radius);

	// particle array and distribution to be used for simulation
	vector<Particle*> parts;
	parts.resize(num_testparts);
	for (int i=0; i<num_testparts; i++)
	{
		parts[i] = new Particle_O();
	}
	Distribution* dist = new Distribution_Hot_O(mars, ref_height, bg_temp);

	// initialize background species
	int num_bgparts = 4;
	Particle* bg_parts[] = {new Particle_O(),
							new Particle_N2(),
							new Particle_CO(),
							new Particle_CO2()};
	double bg_dens[] = {2.64e13, 1.0e13, 3.1e13, 6.68e13};
	double bg_sigs[] = {6.4e-19, 1.85e-18, 1.85e-18, 2.0e-18};
	Distribution_MB* bg_dist = new Distribution_MB(mars, ref_height, bg_temp);
	Background_Species bg_spec(num_bgparts, mars, bg_temp, ref_height, bg_dist, bg_parts, bg_dens, bg_sigs);

	// initialize atmosphere and run simulation
	Atmosphere my_atmosphere(num_testparts, mars, parts, dist, bg_spec, bg_temp, ref_height);
	my_atmosphere.output_velocity_distro(100.0, 150, "/home/rodney/Documents/coronaTest/vdist.out");
	my_atmosphere.run_simulation(dt, timesteps);
	my_atmosphere.output_velocity_distro(100.0, 150, "/home/rodney/Documents/coronaTest/vdist2.out");

*/
	return 0;
}
