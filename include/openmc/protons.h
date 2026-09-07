#ifndef __proton_physics__
#define __proton_physics__

#include <random>
#include <cmath>
//TODO - what RNG to use?? NOT a file local static one pls!!

#include "openmc/nuclide.h"
using Nuclide_t = openmc::Nuclide;
//#include "fake_nuclide.h"

namespace openmc{
static inline std::mt19937 proton_rng {std::random_device {}()};
constexpr double MAX_DEFLECTION = 1.0e-2; // radians

static inline std::uniform_real_distribution<double> uniform_dist {0.0, 1.0};
static inline std::uniform_real_distribution<double> angle_dist {std::cos(MAX_DEFLECTION), 1.0};
static inline std::normal_distribution<double> e_strag(0.0, 1.0);


// IMPORTANT - THIS IS A WIP. A lot of this file is dumb static global state in order to test the MODEL needs before integrating to the codebase proper
//TODO - move some of this into the Nuclide, Material or Particle classes?

inline double mock_random_value(){
    double sample = uniform_dist(proton_rng) * 0.1/0.000668456;
    return sample;
}

inline double random_angle(){
  return angle_dist(proton_rng);
}


/** @brief Inelastic energy loss
 * 
 * Computes the inelastic energy loss per cm using bethe-bloch formula
 * As CURRENTLY implemented this is for a single nuclide in a combined material - the averaged mean-excitation-energy is smuggled in as I and enters non-linearly
 * 
 * @param i_nuclide The index for this nuclide in the global table
 * @param E Initial Energy of the proton in ?????
 * @param I Mean activation energy for current material in ?????
 */
inline double proton_bethe_bloch(int i_nuclide, double E, double I){

    //Access the material base properties from the data table
    //const Nuclide_t& nuclide = *nuclide_list.at(i_nuclide);
    const Nuclide_t& nuclide = *data::nuclides.at(i_nuclide);
    // Betht-bloch contrib for single atom of THIS Nuclide only, summed later
    double mecsq = 0.511;   // mass of electron * speed of light squared, MeV
    double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    E = E / 1e6; // Inside here, expecting MeV
    double betasq = (2 * mpcsq + E) * E / pow(mpcsq + E, 2);
    const double arbitrary_scale = 100;
    return arbitrary_scale*1e6 * 0.3072 * nuclide.Z_ *
             (log(2 * mecsq * betasq / (I * (1 - betasq))) - betasq) /
             (betasq); ///??? In what units??
             ///Removed A from denom as we multiply by this to get mass upstream
    
}

inline double random_straggle(){
  const double arbitrary_scale = 4e17;
  return arbitrary_scale * e_strag(proton_rng);
}
//Random bump - Energy loss or gain due to straggling. Depends on e and material
// and introduced a random Gaussian
inline double energy_straggling_update_sq(double e){
    double alpha = 1 / 137.0;
    double log_hbar = -21 * log(10) + log(4.136) - log(2 * M_PI); // MeV * s
    double log_c = log(29979245800);                              // cm / s
    double log_avogadro = log(6) + 23 * log(10);
    double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    double betasq = (2 * mpcsq + e) * e / pow(mpcsq + e, 2);
    //double log_molecule_density =
        //log(density) + log_avogadro; // molecules / cm^3
    double ret =
        4 * PI * (1 - betasq / 2) / sqrt(1 - betasq) *
        exp(2 * (log(alpha) + log_hbar + log_c));
    return ret;
}
  //Duplicated from SDE code
  inline double energy_straggling_sd(int i_nuclide) {

    //The micro part is just the sum of x * z / a; and can be cached, electrons per average molecule in this material
    //The REST is based on the energy
    const Nuclide_t& nuclide = *data::nuclides.at(i_nuclide);
    return nuclide.Z_ / nuclide.A_;
  }

inline double rutherford_elastic_rate(){
  return 1e4;
}

inline double non_elastic_rate(){
  return 1e4;
}



};
#endif