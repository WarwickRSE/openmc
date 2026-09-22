#include <cmath>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Doing everything in one file for compactness

#include "openmc/proton_cross_sections.h"

//Mock Nuclide to support testing
namespace openmc{
    struct Nuclide{
        double Z_, A_;
        // PROTON TRANSPORT
        CS_1d proton_ne_rate, proton_el_rate; // TODO - share between nuclides??
        CS_2d proton_el_xsec;
        CS_3d proton_ne_xsec;
    };
};

#define __nuclide_included__
#define PI 3.14159265
inline double prn(uint64_t * seed){return 1.0;}
#include "openmc/protons.h"
//------------------------ Basic rate functions ------------------------------------------------

// Check for simple Materials: 1 Water, 2 'Dense' carbon, 3 'bone'

TEST_CASE("Proton Rates in Single Nuclide"){

    //Fake Nuclide
    openmc::Nuclide N;
    openmc::read_proton_data(N, "/media/raid/MathRadData/protons/", "H1");
    
    N.Z_ = 1;
    N.A_ = 1;

    std::vector<double> energies{200, 150, 100, 75, 40};
    for(int i = 0; i < 5; i++){
      auto E = energies[i]*1e6;
      auto I = 70.0;
      std::cerr<<openmc::proton_sde::proton_bethe_bloch(N, E, I)<<std::endl;
    }

    // Checking rate at the listed energies

}