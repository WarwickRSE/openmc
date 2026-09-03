#ifndef __proton_physics__
#define __proton_physics__

#include "nuclide.h"
#include <random>
#include <cmath>
//TODO - what RNG to use?? NOT a file local static one pls!!

namespace openmc{
static inline std::mt19937 proton_rng {std::random_device {}()};
static inline std::uniform_real_distribution<double> uniform_dist {0.0, 1.0};

double mock_x_section(){
    return 0.1/0.000668456; // A value which makes something happen
}

//TODO - move some of this into the Nuclide class?

double mock_random_value(){
    double sample = uniform_dist(proton_rng) * 0.1/0.000668456;
    return sample;
}


double proton_bethe_bloch(int i_nuclide, double E, double I){

    //Access the material base properties from the data table
    const Nuclide& nuclide = *data::nuclides.at(i_nuclide);
    // Betht-bloch contrib for THIS Nuclide only, summed later
    double mecsq = 0.511;   // mass of electron * speed of light squared, MeV
    double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    double betasq = (2 * mpcsq + E) * E / pow(mpcsq + E, 2);
    
    //TODO - density* mass_fraction == atom_density??
    return 0.3072 * nuclide.Z_ *
             (log(2 * mecsq * betasq / (I * (1 - betasq))) - betasq) /
             (betasq * nuclide.A_); // MeV / cm
    
}

};
#endif