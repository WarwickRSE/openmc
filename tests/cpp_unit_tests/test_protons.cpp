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

constexpr double eps_calc = 1e-5;

// Check for simple Materials: 1 Water, 2 'Dense' carbon, 3 'bone'

TEST_CASE("Proton Rates in Single Nuclide"){
    //Checking partial contributions from selected Nuclides

    //Fake Nuclide - Hydrogen
    openmc::Nuclide H1;
    openmc::read_proton_data(&H1.proton_el_rate, &H1.proton_ne_rate, &H1.proton_el_xsec, &H1.proton_ne_xsec, "/media/raid/MathRadData/protons/", "H1");
    H1.Z_ = 1;
    H1.A_ = 1.008;
    
    //Fake Nuclide - Oxygen
    openmc::Nuclide O16;
    openmc::read_proton_data(&O16.proton_el_rate, &O16.proton_ne_rate, &O16.proton_el_xsec, &O16.proton_ne_xsec, "/media/raid/MathRadData/protons/", "O16");
    O16.Z_ = 8;
    O16.A_ = 15.999;

    //Fake Nuclide - Carbon
    openmc::Nuclide C12;
    openmc::read_proton_data(&C12.proton_el_rate, &C12.proton_ne_rate, &C12.proton_el_xsec, &C12.proton_ne_xsec, "/media/raid/MathRadData/protons/", "C12");
    C12.Z_ = 6;
    C12.A_ = 12.011;


    std::vector<double> energies{200, 150, 100, 75, 40};

    // Created by test code
    std::vector<double> ref_loss_H1{8.09683, 9.81462, 13.1409, 16.3419, 26.8417};
    for(int i = 0; i < 5; i++){
      auto E = energies[i]*1e6;
      auto I = 75.0; //Correct for e.g. H in Water
      auto loss = openmc::proton_sde::proton_bethe_bloch(H1, E, I)/1e6;
      REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_loss_H1[i], eps_calc));
    }

    // Created by test code
    std::vector<double> ref_loss_O16{4.04867, 4.90762, 6.57085, 8.17145, 13.4217};
    for(int i = 0; i < 5; i++){
      auto E = energies[i]*1e6;
      auto I = 75.0; //Correct for e.g. H in Water
      // Divide by A_ as we canclled this factor
      auto loss = openmc::proton_sde::proton_bethe_bloch(O16, E, I)/1e6/O16.A_;
      REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_loss_O16[i], eps_calc));
    }

    std::vector<double> ref_loss_C12{3.98478, 4.82798, 6.45965, 8.02871, 13.1688};
    for(int i = 0; i < 5; i++){
      auto E = energies[i]*1e6;
      auto I = 85.0; //Correct for e.g. H in Water
      // Divide by A_ as we canclled this factor
      auto loss = openmc::proton_sde::proton_bethe_bloch(C12, E, I)/1e6/C12.A_;
      REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_loss_C12[i], eps_calc));
    }


}