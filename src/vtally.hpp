/*
 * vtally.hpp
 *
 *  Created on: Jun 19, 2020
 *      Author: rodney
 */

#ifndef VTALLY_HPP_
#define VTALLY_HPP_

#include "Atmosphere.hpp"
#include "Particle.hpp"

class Vtally
{

private:
  static const int num_vel_bins = 90; // set number of velocity bins here
  static const double LOS_velocity_spare;
  static const int vel_bin;
  static constexpr double vel_max = 45.0;
  static constexpr double vel_min = 0.0;
  const double dv = (vel_max-vel_min)/num_vel_bins;
  static constexpr double dx = 10; // if particle's x coordinate is between chosen altitude and chosen altitude + dx, add it to the count (km)
  static const int no_LOS_angle_bins = 12; // no of LOS angles over which vel dist is averaged
  const double angle_spacing = constants::pi/no_LOS_angle_bins;
  static constexpr double w = 10e5; // width of line of sight, in y-z plane (cm)
  vector<double> chosen_alts = {150, 300, 1000}; //altitudes (km) at which we simulate a LOS
  double weighting_factor;
  double radius_km;
  int bin;
  double vel_LOS;
  
  vector<vector<double>> vtally_matrix; // rows = no. of velocity bins; columns = no. of altitudes (from corona3d.cfg); contains counts and is added to at each timestep.
  vector<vector<double>> cos_sin_tan_angle; // rows: [0] = cos(angle); [1] = sin(angle); [2] = tan(angle); columns: one for each angle


public:
  Vtally(const int &num_EDFs, const vector<int> &EDF_alts, const double &vtally_x, const double &vtally_dx,  const double &vtally_w, const double &dt, const double &rate, const int &num_parts, const double &radius); //, double rate);

  void update_vtally(const shared_ptr<Particle> &p);
  void calculate_LOS_velocity(const shared_ptr<Particle> &p, const int &j);
  int choose_vel_bin(const double &vel) const;
  bool is_inside(const shared_ptr<Particle> &p, const double &alt) const;
  void record_vtallies(const string &output_dir);
};

#endif /* VTALLY_HPP_ */
