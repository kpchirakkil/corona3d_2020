/*
 * vtally.cpp
 *
 *  Created on: Jun 12, 2023
 *      Author: bethan
 */
#include "Atmosphere.hpp"
#include "Vtally.hpp"
#include <iostream>

Vtally::Vtally(const int &num_EDFs, const vector<int> &EDF_alts, const double &vtally_x, const double &vtally_dx,  const double &vtally_w, const double &dt, const double &rate, const int &num_parts, const double &radius){ //
  radius_km = radius/1e5;
  weighting_factor = dt*rate/(no_LOS_angle_bins*num_parts*(dx*1e5)*w);

  vtally_matrix.resize(num_vel_bins);
  for (int i=0; i<num_vel_bins; i++)
    {
      vtally_matrix[i].resize(num_EDFs);
      for (int j=0; j<num_EDFs; j++)
      {
        vtally_matrix[i][j] = 0.0;
      }
    }
  cos_sin_tan_angle.resize(3); // array for cos(angle), sin(angle), and tan(angle) values
  for (int i=0; i<3; i++){
    cos_sin_tan_angle[i].resize(no_LOS_angle_bins);
  }
  for (int j=0; j<no_LOS_angle_bins; j++){
    cos_sin_tan_angle[0][j] = cos(j*angle_spacing);
    cos_sin_tan_angle[1][j] = sin(j*angle_spacing);
    cos_sin_tan_angle[2][j] = tan(j*angle_spacing);
    
   }
}

void Vtally::update_vtally(const shared_ptr<Particle> &p){

  for (int j=0; j< chosen_alts.size(); j++)
    {
      if (is_inside(p, chosen_alts[j]))
	{
	  calculate_LOS_velocity(p, j);
	}
    }
}

void Vtally::calculate_LOS_velocity(const shared_ptr<Particle> &p, const int &j){
  const double y = p->get_y(); //(cm)
  const double z = p->get_z(); //(cm)
  const double vy = p->get_vy(); //(cm)
  const double vz = p->get_vz(); //(cm)
  
  for (int l=0; l<no_LOS_angle_bins; l++){ // make this more efficient - probably don't need this loop... maybe? Probably a more efficient way
    const double angle = l*angle_spacing;
	 if (angle == 0.0){
	   if (z >= -1*w/2 and z<= w/2){
	     vtally_matrix[choose_vel_bin(abs(vy))][j] += weighting_factor;
	     
	   }
	 }
	 else if (angle == constants::pi/2){
	   if (y>= -1*w/2 and y<= w/2){
	     vtally_matrix[choose_vel_bin(abs(vz))][j] += weighting_factor;
	     
	   }
         }
	 else {
	   if (y >= (-1*z)/cos_sin_tan_angle[2][l] - w/(2*cos_sin_tan_angle[1][l]) && y <= (-1*z)/cos_sin_tan_angle[2][l] + w/(2*cos_sin_tan_angle[1][l])){
	     vel_LOS = abs(vy*cos_sin_tan_angle[0][l] - vz*cos_sin_tan_angle[1][l]); // this is equal to the component of vy along the LOS + the component of vz along the LOS = vy*cos(angle) = vz*sin(angle)
	     vtally_matrix[choose_vel_bin(vel_LOS)][j] += weighting_factor;
	 }
	 }
	 
  }
  }

int Vtally::choose_vel_bin(const double &LOS_velocity) const {
  return std::floor(((LOS_velocity*1e-5)-vel_min)/((vel_max-vel_min)/num_vel_bins));
}


// find out whether particle i has x co-ordinate at chosen altitude (between chosen altitude and 1 km greater)
bool Vtally::is_inside(const shared_ptr<Particle> &p, const double &alt) const {  //(int i)
 const double alt_p = p->get_x()*1e-5 - radius_km; // in km
 return (alt_p >= alt && alt_p < alt+dx);
}

void Vtally::record_vtallies(const string &output_dir){
  ofstream vtally_out;

  for (std::size_t k = 0; k < chosen_alts.size(); k++){
    int int_alt = chosen_alts[k];
    vtally_out.open(output_dir + "vtally_" + to_string(int_alt) + "km.csv");
    vtally_out << "#velocity_bin_no,min_velocity_in_bin_[km/s],column_density_[cm-2]\n";
    // Copy counts for positive velocities (that were tallied) to negative velocities, as it is a symmetrical distribution.
    for (int j=0; j < num_vel_bins; j++){
      vtally_out << num_vel_bins - 1 - j << "," << (-1*vel_max) + j*((vel_max-vel_min)/num_vel_bins) << "," << vtally_matrix[num_vel_bins - 1 - j][k] << "\n";
    }
    // Record counts for positive velocities
    for (int j=0; j < num_vel_bins; j++){
      vtally_out << j << "," << vel_min + j*((vel_max-vel_min)/num_vel_bins) << "," << vtally_matrix[j][k] <<"\n"; // here, j is the vel_bin count and i is the altitude count
  }
    vtally_out.close();
  }
}

  
