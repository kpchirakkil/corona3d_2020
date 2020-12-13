/*
 * Atmosphere.hpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#ifndef ATMOSPHERE_HPP_
#define ATMOSPHERE_HPP_

#include <vector>
#include <iomanip>
#include "Background_Species.hpp"
#include "Distribution_Hot_H.hpp"
#include "Distribution_Hot_O.hpp"
#include "Distribution_Import.hpp"
#include "Distribution_MB.hpp"
#include "Common_Functions.hpp"
using namespace std;

class Atmosphere {
public:
	Atmosphere(int n, int num_to_trace, string trace_output_dir, Planet p, vector<shared_ptr<Particle>> parts, shared_ptr<Distribution> dist, Background_Species bg, int num_EDFs, int EDF_alts[]);
	virtual ~Atmosphere();

	void output_positions(string datapath);
	void output_altitude_distro(double bin_width, string datapath);
	void output_velocity_distro(double bin_width, string datapath);
	void run_simulation(double dt, int num_steps, double lower_bound, double upper_bound, double avg_thermal_v, int print_status_freq, int output_pos_freq, string output_pos_dir, string output_stats_dir);

private:
	int num_parts;                      // number of particles initially spawned
	int num_traced;                     // number of particles to output trace data on
	string trace_dir;                   // directory to output particle trace data to
	int active_parts;                   // number of active particles
	Planet my_planet;                   // contains planet mass and radius
	vector<shared_ptr<Particle>> my_parts;         // array of particles to be tracked
	shared_ptr<Distribution> my_dist;              // distribution class to initialize particles
	Background_Species bg_species;      // background species used for collisions
	vector<int> traced_parts;           // indices of randomly selected trace particles

	vector<int> stats_dens_counts;
	int stats_num_EDFs;
	vector<int> stats_EDF_alts;
	vector<vector<int>> stats_EDFs;
	void update_stats();
	void output_stats(double dt, double rate, int total_parts, string output_dir);

	// output test particle trace data for selected particles
	void output_collision_data();
	void output_trace_data();
};

#endif /* ATMOSPHERE_HPP_ */
