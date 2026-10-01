#ifndef __cross_sections__
#define __cross_sections__

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace openmc{

class Nuclide;

constexpr inline double divide_by_zero_tol = 1e-7; //Allowable minimum denominator to avoid divide-by-zero

//! Convert from centre-of-mass frame to lab frame
//! Assumes 'other' nuclide is Hydrogen (for the mass)
inline constexpr double hydrogen_cm_to_lab(const double ang, const double E) {
    double mp = 938.346;
    double E2mp = E + 2.0 * mp;
    double p = sqrt(E * E2mp);
    double u = p / E2mp;
    double g = 1.0 / sqrt(1.0 - u * u);
    double e_pr = E + mp;
    double v_ratio = u * (e_pr - u * p) / (p - u * e_pr);
    return atan2(sin(ang), (g * (v_ratio - cos(ang))));
}

//! Linear interpolation helper
//! right_offset is the distance between x_left and target, i..e the distance into the target cell
//! NOTE: if dx is too small, this returns right_val
inline constexpr double interp(const double left_val, const double right_val, const double left_offset, const double dx){

  if(dx > divide_by_zero_tol){
    const double frac = left_offset/dx;
    return (left_val * (1.0 - frac) + right_val * frac);
  }else{
    return right_val;
  }
}

struct CS_1d {

  CS_1d(const std::string filename) : energy(), rate() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    std::stringstream iss;

    getline(file, line);
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }

    getline(file, line);
    std::stringstream iss2;
    iss2 << line;
    while (getline(iss2, token, ' ')) {
      rate.push_back(atof(token.c_str()));
    }
    file.close();
    if(energy.size() != rate.size()){
      throw std::runtime_error("Data file "+filename+" invalid. Energies and rates do not match in length");
    }
  }
  CS_1d(const std::string filename, const double cuttoff) : energy(), rate() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    double tmp_val, tmp_val_old = 0, lin_inter_val = 0;
    int tmp_count, tmp_count_2;
    bool lin_inter_bool;
    while (getline(file, line)) {
      tmp_count = 0;
      lin_inter_bool = true;
      std::stringstream iss2;
      iss2 << line;
      while (getline(iss2, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_val > cuttoff) {
          tmp_count++;
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          lin_inter_val = (cuttoff - tmp_val_old) / (tmp_val - tmp_val_old);
          lin_inter_bool = false;
        }
      }
      tmp_count_2 = 0;
      lin_inter_bool = true;
      getline(file, line);
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_count_2 < tmp_count) {
          tmp_count_2++;
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          rate.push_back(tmp_val * lin_inter_val +
                         (1 - lin_inter_val) * tmp_val_old);
          lin_inter_bool = false;
        }
      }
      if (lin_inter_bool) {
        rate.push_back(tmp_val_old);
      }
    }
    file.close();
    if(energy.size() != rate.size()){
      throw std::runtime_error("Data file "+filename+" invalid. Energies and rates do not match in length");
    }
  }

  CS_1d(const CS_1d &other) : energy(other.energy), rate(other.rate) {}
  CS_1d & operator=(const CS_1d & other){energy=other.energy;rate=other.rate; return *this;}

  CS_1d(const std::string filename, const double cuttoff,
        const double back_cuttoff)
      : energy(), rate() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    double lab_ang_cutoff, tmp_val, tmp_val_old = 0, top_rate, bottom_rate,
                                    total_rate, lin_inter_val = 0;
    int tmp_count, tmp_count_2, tmp_count_back, tmp_count_back_2,
        energy_index = 0;
    bool lin_inter_bool;
    while (getline(file, line)) {
      lab_ang_cutoff = hydrogen_cm_to_lab(back_cuttoff, energy[energy_index]);
      energy_index++;
      tmp_count = 0;
      tmp_count_back = 0;
      lin_inter_bool = true;
      std::stringstream iss2;
      iss2 << line;
      while (getline(iss2, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_val > lab_ang_cutoff) {
          tmp_count++;
          tmp_count_back++;
        } else if (tmp_val > cuttoff) {
          tmp_count++;
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          lin_inter_val = (cuttoff - tmp_val_old) / (tmp_val - tmp_val_old);
          lin_inter_bool = false;
        }
      }
      tmp_count_2 = 0;
      tmp_count_back_2 = 0;
      lin_inter_bool = true;
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_count_back_2 < tmp_count_back) {
          tmp_count_back_2++;
          tmp_count_2++;
        } else if (tmp_count_2 < tmp_count) {
          tmp_count_2++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          tmp_vec.push_back(tmp_val * lin_inter_val +
                            (1 - lin_inter_val) * tmp_val_old);
          lin_inter_bool = false;
        }
      }
      top_rate = tmp_vec.back();
      bottom_rate = tmp_vec.front();
      total_rate = top_rate - bottom_rate;
      rate.push_back(total_rate);
    }
    file.close();
    if(energy.size() != rate.size()){
      throw std::runtime_error("Data file "+filename+" invalid. Energies and rates do not match in length");
    }
  }

  CS_1d() : energy(), rate() {}

  double evaluate(const double e) const {
    if (energy.size() > 0) {
      // Find cell
      int i_r = std::distance(energy.begin(), std::lower_bound(energy.begin(), energy.end(), e));
      if (i_r == 0){
        return rate[0];
      } else if (i_r  == int(rate.size())){
        return rate.back();
      } else {
       // Interpolate within cell
        return interp(rate[i_r-1], rate[i_r], (e - energy[i_r-1]), (energy[i_r] - energy[i_r-1]));
      }
    }else{
      return 0.0;
    }
  }
  std::vector<double> energy, rate;
};

struct CS_3d {

  CS_3d(const std::string filename) : energy(), exit_energy(), cdf(), rvalue() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    while (getline(file, line)) {
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      exit_energy.push_back(tmp_vec);
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss4;
      iss4 << line;
      while (getline(iss4, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      cdf.push_back(tmp_vec);
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss5;
      iss5 << line;
      while (getline(iss5, token, ' ')) {
        tmp_vec.push_back(atof(token.c_str()));
      }
      rvalue.push_back(tmp_vec);
    }
    file.close();
    if (energy.size() != exit_energy.size() || energy.size() != cdf.size() ||
        energy.size() != rvalue.size()) {
      throw std::runtime_error("Data file " + filename +
                               " invalid. Energy, exit-energy, CDF, and "
                               "r-value slice counts do not match");
    }
    for (std::size_t i = 0; i < energy.size(); ++i) {
      if (exit_energy[i].size() != cdf[i].size() ||
          exit_energy[i].size() != rvalue[i].size()) {
        throw std::runtime_error("Data file " + filename + " invalid. "
                                 "Energy slice " + std::to_string(i) +
                                 " has mismatched exit-energy, CDF, and "
                                 "r-value lengths");
      }
    }
  }

  CS_3d(const CS_3d &other)
      : energy(other.energy), exit_energy(other.exit_energy), cdf(other.cdf),
        rvalue(other.rvalue) {}

  CS_3d() : energy(), exit_energy(), cdf(), rvalue() {}
  
  CS_3d & operator=(const CS_3d & other){energy=other.energy; exit_energy=other.exit_energy; cdf=other.cdf;rvalue=other.rvalue; return *this;}


  //! Sample from selected slice
  // ! Interpolate energy and r based on position in cdf
  std::pair<double, double> sample_from_vector(const int energy_index, const double u) const {

    const auto & cdf_slice = cdf[energy_index];
    const auto & e_slice = exit_energy[energy_index];
    const auto & r_slice = rvalue[energy_index];
    const int i_d = std::distance(cdf_slice.begin(), std::lower_bound(cdf_slice.begin(), cdf_slice.end(), u));

    if (i_d == 0) {
      return {e_slice[0], r_slice[0]};
    } else if (i_d == int(cdf_slice.size())) {
      return {e_slice.back(), r_slice.back()};
    } else {
      const auto inc = u - cdf_slice[i_d-1];
      const auto dx = cdf_slice[i_d] - cdf_slice[i_d-1];
      auto e = interp(e_slice[i_d-1], e_slice[i_d], inc, dx);
      auto r = interp(r_slice[i_d-1], r_slice[i_d], inc, dx);
      return {e, r};
    }
  }

  std::pair<double, double> sample(const double e, double u) const {
    const int i_e = std::distance(energy.begin(), std::lower_bound(energy.begin(), energy.end(), e));
    if (i_e == 0) {
      return sample_from_vector(0, u);
    } else if (i_e == int(energy.size())) {
      return sample_from_vector(i_e - 1, u);
    } else {
      // NOTE: in case the energy axis does not meet tolerance, this will
      // do one needless interpolation. However, that should be a rare case
      auto left_sample  = sample_from_vector(i_e - 1, u);
      auto right_sample = sample_from_vector(i_e,     u);
      auto dx = energy[i_e] - energy[i_e - 1];
      auto diff = e - energy[i_e - 1];
      auto e = interp(left_sample.first , right_sample.first , diff, dx);
      auto r = interp(left_sample.second, right_sample.second, diff, dx);
      return {e, r};
    }
  }

  std::vector<double> energy;
  std::vector<std::vector<double>> exit_energy, cdf, rvalue;
};

struct CS_2d {

  CS_2d(const std::string filename, const double cuttoff)
      : energy(), exit_angle(), cdf() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    double tmp_val, tmp_val_old = 0, total_rate, lin_inter_val = 0;
    int tmp_count, tmp_count_2;
    bool lin_inter_bool;
    while (getline(file, line)) {
      tmp_vec.clear();
      tmp_count = 0;
      lin_inter_bool = true;
      std::stringstream iss2;
      iss2 << line;
      while (getline(iss2, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_val > cuttoff) {
          tmp_count++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          lin_inter_val = (cuttoff - tmp_val_old) / (tmp_val - tmp_val_old);
          lin_inter_bool = false;
          tmp_vec.push_back(cuttoff);
        }
      }
      exit_angle.push_back(tmp_vec);
      tmp_count_2 = 0;
      lin_inter_bool = true;
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_count_2 < tmp_count) {
          tmp_count_2++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          tmp_vec.push_back(tmp_val * lin_inter_val +
                            (1 - lin_inter_val) * tmp_val_old);
          lin_inter_bool = false;
        }
      }
      total_rate = tmp_vec.back();
      for (double &i : tmp_vec) {
        i /= total_rate;
      }
      cdf.push_back(tmp_vec);
    }
    file.close();
    if (energy.size() != exit_angle.size() || energy.size() != cdf.size()) {
      throw std::runtime_error("Data file " + filename +
                               " invalid. Energy, exit-angle, and CDF "
                               "slice counts do not match");
    }
    for (std::size_t i = 0; i < energy.size(); ++i) {
      if (exit_angle[i].size() != cdf[i].size()) {
        throw std::runtime_error("Data file " + filename + " invalid. "
                                 "Energy slice " + std::to_string(i) +
                                 " has mismatched exit-angle and CDF lengths");
      }
    }
  }

  CS_2d(const std::string filename, const double cuttoff,
        const double back_cuttoff)
      : energy(), exit_angle(), cdf() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    getline(file, line);
    std::stringstream iss;
    iss << line;
    while (getline(iss, token, ' ')) {
      energy.push_back(atof(token.c_str()));
    }
    std::vector<double> tmp_vec;
    double lab_ang_cutoff, tmp_val, tmp_val_old = 0, top_rate, bottom_rate,
                                    total_rate, lin_inter_val = 0;
    int tmp_count, tmp_count_2, tmp_count_back, tmp_count_back_2,
        energy_index = 0;
    bool lin_inter_bool;
    while (getline(file, line)) {
      lab_ang_cutoff = hydrogen_cm_to_lab(back_cuttoff, energy[energy_index]);
      energy_index++;
      tmp_vec.clear();
      tmp_count = 0;
      tmp_count_back = 0;
      lin_inter_bool = true;
      std::stringstream iss2;
      iss2 << line;
      while (getline(iss2, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_val > lab_ang_cutoff) {
          tmp_count++;
          tmp_count_back++;
        } else if (tmp_val > cuttoff) {
          tmp_count++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          lin_inter_val = (cuttoff - tmp_val_old) / (tmp_val - tmp_val_old);
          lin_inter_bool = false;
          tmp_vec.push_back(cuttoff);
        }
      }
      exit_angle.push_back(tmp_vec);
      tmp_count_2 = 0;
      tmp_count_back_2 = 0;
      lin_inter_bool = true;
      getline(file, line);
      tmp_vec.clear();
      std::stringstream iss3;
      iss3 << line;
      while (getline(iss3, token, ' ')) {
        tmp_val = atof(token.c_str());
        if (tmp_count_back_2 < tmp_count_back) {
          tmp_count_back_2++;
          tmp_count_2++;
        } else if (tmp_count_2 < tmp_count) {
          tmp_count_2++;
          tmp_vec.push_back(tmp_val);
          tmp_val_old = tmp_val;
        } else if (lin_inter_bool) {
          tmp_vec.push_back(tmp_val * lin_inter_val +
                            (1 - lin_inter_val) * tmp_val_old);
          lin_inter_bool = false;
        }
      }
      top_rate = tmp_vec.back();
      bottom_rate = tmp_vec.front();
      total_rate = top_rate - bottom_rate;
      for (double &i : tmp_vec) {
        i = (i - bottom_rate) / total_rate;
      }
      cdf.push_back(tmp_vec);
    }
    file.close();
    if (energy.size() != exit_angle.size() || energy.size() != cdf.size()) {
      throw std::runtime_error("Data file " + filename +
                               " invalid. Energy, exit-angle, and CDF "
                               "slice counts do not match");
    }
    for (std::size_t i = 0; i < energy.size(); ++i) {
      if (exit_angle[i].size() != cdf[i].size()) {
        throw std::runtime_error("Data file " + filename + " invalid. "
                                 "Energy slice " + std::to_string(i) +
                                 " has mismatched exit-angle and CDF lengths");
      }
    }
  }

  CS_2d(const CS_2d &other)
      : energy(other.energy), exit_angle(other.exit_angle), cdf(other.cdf) {}

  CS_2d() : energy(), exit_angle(), cdf() {}
  CS_2d & operator=(const CS_2d & other){energy=other.energy; exit_angle=other.exit_angle; cdf=other.cdf; return *this;} 

  //! Sample from selected slice
  // ! Interpolate energy and r based on position in cdf
  double sample_from_vector(const int energy_index, const double u) const {

    const auto & cdf_slice = cdf[energy_index];
    const auto & a_slice = exit_angle[energy_index];
    const int i_d = std::distance(cdf_slice.begin(), std::lower_bound(cdf_slice.begin(), cdf_slice.end(), u));

    if (i_d == 0) {
      return a_slice[0];
    } else if (i_d == int(cdf_slice.size())) {
      return a_slice.back();
    } else {
      const auto inc = u - cdf_slice[i_d-1];
      const auto dx = cdf_slice[i_d] - cdf_slice[i_d-1];
      return interp(a_slice[i_d-1], a_slice[i_d], inc, dx);
    }
  }
  
  double sample(const double e, double u) const {
    const int i_e = std::distance(energy.begin(), std::lower_bound(energy.begin(), energy.end(), e));
    if (i_e == 0) {
      return sample_from_vector(0, u);
    } else if (i_e == int(energy.size())) {
      return sample_from_vector(i_e - 1, u);
    } else {
      // NOTE: in case the energy axis does not meet tolerance, this will
      // do one needless interpolation. However, that should be a rare case
      auto left_sample  = sample_from_vector(i_e - 1, u);
      auto right_sample = sample_from_vector(i_e,     u);
      auto dx = energy[i_e] - energy[i_e - 1];
      auto diff = e - energy[i_e - 1];
      return interp(left_sample, right_sample, diff, dx);
    }
  }

  std::vector<double> energy;
  std::vector<std::vector<double>> exit_angle, cdf;
};

// Reading functions
//! Convert from an Isotope name to the long-form element name
std::string letter_to_string(std::string sym);
//! Read the cross sections for specified name, and data file path
void read_proton_data(CS_1d * ne_rate, CS_1d* el_rate, CS_2d * el_xsec, CS_3d * ne_xsec, std::string path, std::string name);

}
#endif
