#include <cmath>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Doing everything in one file for compactness

#include "openmc/constants.h"
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

// Threshold for equality, keeping in mind that we might be using
// slightly different calculation breakdowns so do NOT expect full
// FP equality
constexpr double eps_calc = 1e-5;

//TODO - get from somewhere... like char* proton_data_path = std::getenv("OPENMC_PROTON_DATA");
const std::string data_path = "/media/raid/MathRadData/protons/";

TEST_CASE("Proton Energy Rates in Single Nuclide"){
    //Checking partial contributions from selected Nuclides

    //Fake Nuclide - Hydrogen
    openmc::Nuclide H1;
    H1.Z_ = 1;
    H1.A_ = 1.008;
    
    //Fake Nuclide - Oxygen
    openmc::Nuclide O16;
    O16.Z_ = 8;
    O16.A_ = 15.999;

    //Fake Nuclide - Carbon
    openmc::Nuclide C12;
    C12.Z_ = 6;
    C12.A_ = 12.011;


    std::vector<double> energies{200, 150, 100, 75, 40};

    SECTION("Hydrogen Bethe Bloch"){
        // Created by test code
      std::vector<double> ref_loss_H1{8.09683, 9.81462, 13.1409, 16.3419, 26.8417};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        auto I = 75.0; //Correct for e.g. H in Water
        auto loss = openmc::proton_sde::proton_bethe_bloch(H1, E, I)/1e6;
        REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_loss_H1[i], eps_calc));
       }
    }
    SECTION("Oxygen Bethe Bloch"){
      // Created by test code
      std::vector<double> ref_loss_O16{4.04867, 4.90762, 6.57085, 8.17145, 13.4217};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        auto I = 75.0; //Correct for e.g. H in Water
        // Divide by A_ as we canclled this factor
        auto loss = openmc::proton_sde::proton_bethe_bloch(O16, E, I)/1e6/O16.A_;
        REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_loss_O16[i], eps_calc));
      }
    }
    SECTION("Carbon Bethe Bloch"){
      std::vector<double> ref_loss_C12{3.98478, 4.82798, 6.45965, 8.02871, 13.1688};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        auto I = 85.0; //Made up, matched to reference generator
        // Divide by A_ as we canclled this factor
        auto loss = openmc::proton_sde::proton_bethe_bloch(C12, E, I)/1e6/C12.A_;
        REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_loss_C12[i], eps_calc));
      }
    }


    //NOTE: for straggling we do partially reproduce the logic from physics.cpp, particle.cpp
    // and material.cpp AND assume a density of 1.0
    // HOWEVER we will check everything again later. Here we just accept it, because it
    // avoids having to set up so much infrastructure, and stub out the RNG
    //Consider this test to be a cross-reference for what the core should do and a check on the basic
    // Nuclide dependency. Note again we cancel a factor A_
    SECTION("Hydrogen Energy Straggling"){
        // Created by test code
      std::vector<double> ref_strag_H1{0.43796, 0.426611, 0.415482, 0.410006, 0.402444};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        double distance = 1.0;
        auto loss = std::sqrt(openmc::proton_sde::energy_straggling_sd(H1) * exp(openmc::proton_sde::log_avogadro) / H1.A_ * openmc::proton_sde::energy_straggling_update_sq(E) * distance);
        REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_strag_H1[i], eps_calc));
       }
    }
    SECTION("Oxygen Energy Straggling"){
        // Created by test code
      std::vector<double> ref_strag_O16{0.310931, 0.302873, 0.294972, 0.291084, 0.285716};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        double distance = 1.0;
        auto loss = std::sqrt(openmc::proton_sde::energy_straggling_sd(O16) * exp(openmc::proton_sde::log_avogadro) /O16.A_ * openmc::proton_sde::energy_straggling_update_sq(E) * distance);
        REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_strag_O16[i], eps_calc));
       }
    }
    SECTION("Carbon Energy Straggling"){
        // Created by test code
      std::vector<double> ref_strag_C12{0.310779, 0.302725, 0.294828, 0.290942, 0.285576};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        double distance = 1.0;
        auto loss = std::sqrt(openmc::proton_sde::energy_straggling_sd(C12) * exp(openmc::proton_sde::log_avogadro) /C12.A_ * openmc::proton_sde::energy_straggling_update_sq(E) * distance);
        REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_strag_C12[i], eps_calc));
       }
    }

}

TEST_CASE("Proton Cross Sections from File"){
    //Checking partial contributions from selected Nuclides
    // We have to reconstruct the density component and I think we have an N_AVOGADRO stray in these

    //Fake Nuclide - Hydrogen
    openmc::Nuclide H1;
    openmc::read_proton_data(&H1.proton_el_rate, &H1.proton_ne_rate, &H1.proton_el_xsec, &H1.proton_ne_xsec, data_path, "H1");
    H1.Z_ = 1;
    H1.A_ = 1.008;
    
    //Fake Nuclide - Oxygen
    openmc::Nuclide O16;
    openmc::read_proton_data(&O16.proton_el_rate, &O16.proton_ne_rate, &O16.proton_el_xsec, &O16.proton_ne_xsec, data_path, "O16");
    O16.Z_ = 8;
    O16.A_ = 15.999;

    //Fake Nuclide - Carbon
    openmc::Nuclide C12;
    openmc::read_proton_data(&C12.proton_el_rate, &C12.proton_ne_rate, &C12.proton_el_xsec, &C12.proton_ne_xsec, data_path, "C12");
    C12.Z_ = 6;
    C12.A_ = 12.011;

    std::vector<double> energies{200, 150, 100, 75, 40};

    SECTION("Hydrogen Elastic Rate"){
        // Created by test code
      std::vector<double> ref_el_H1{0.113483, 0.113483, 0.157278, 0.171171, 0.427228};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        double atom_density = 1.0 * openmc::N_AVOGADRO / H1.A_;
        auto rate = openmc::proton_sde::rutherford_elastic_rate(H1, E) * atom_density;
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_el_H1[i], eps_calc));
       }
    }
    SECTION("Oxygen Elastic Rate"){
        // Created by test code
      std::vector<double> ref_el_O16{0.00693154, 0.00693154, 0.0150371, 0.0255221, 0.0752361};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        double atom_density = 1.0 * openmc::N_AVOGADRO / O16.A_;
        auto rate = openmc::proton_sde::rutherford_elastic_rate(O16, E) * atom_density;
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_el_O16[i], eps_calc));
       }
    }
    SECTION("Carbon Elastic Rate"){
        // Created by test code
      std::vector<double> ref_el_C12{0.0134684, 0.0134684, 0.0287689, 0.0477845, 0.128913};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        double atom_density = 2.0 * openmc::N_AVOGADRO / C12.A_;
        auto rate = openmc::proton_sde::rutherford_elastic_rate(C12, E) * atom_density;
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_el_C12[i], eps_calc));
       }
    }

    SECTION("Hydrogen Non-Elastic Rate"){
      // Value 0, but expect no-throwing
      REQUIRE_NOTHROW(openmc::proton_sde::non_elastic_rate(H1, 0.0));
      REQUIRE(openmc::proton_sde::non_elastic_rate(H1, 1.0e6) == 0.0);
    }
    SECTION("Oxygen Non-Elastic Rate"){
        // Created by test code
      std::vector<double> ref_ne_O16{0.0110632, 0.0110632, 0.0111382, 0.0121133, 0.0168386};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        double atom_density = 1.0 * openmc::N_AVOGADRO / O16.A_;
        auto rate = openmc::proton_sde::non_elastic_rate(O16, E) * atom_density;
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_ne_O16[i], eps_calc));
       }
    }
    SECTION("Carbon Non-Elastic Rate"){
        // Created by test code
      std::vector<double> ref_ne_C12{0.0222146, 0.0222146, 0.0226792, 0.0268254, 0.0362666};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        double atom_density = 2.0 * openmc::N_AVOGADRO / C12.A_;
        auto rate = openmc::proton_sde::non_elastic_rate(C12, E) * atom_density;
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_ne_C12[i], eps_calc));
       }
    }
 
  }

  TEST_CASE("Proton Moliere Scattering"){
    //Checking partial contributions from selected Nuclides
    //Once again we are reproducing some of the transforms, but this helps cross-check the core
    // calculation
    //Fake Nuclide - Hydrogen
    openmc::Nuclide H1;
    H1.Z_ = 1;
    H1.A_ = 1.008;
    
    //Fake Nuclide - Oxygen
    openmc::Nuclide O16;
    O16.Z_ = 8;
    O16.A_ = 15.999;

    //Fake Nuclide - Carbon
    openmc::Nuclide C12;
    C12.Z_ = 6;
    C12.A_ = 12.011;

    std::vector<double> energies{200, 150, 100, 75, 40};

    const double fixed_step = 0.05;
    SECTION("Hydrogen Small Angle rate"){

      std::vector<double> ref_sa_H1{0.00394011, 0.00519039, 0.00771194, 0.0102573, 0.019297};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        double density = 1.0;
        auto tmp = openmc::proton_sde::moliere_scattering_precomp(H1, E);
        //Single nuclide, so no need to sum anything
        auto sd = openmc::proton_sde::moliere_transform(E, tmp.first, tmp.second, density);
        auto rate =  std::sqrt((1.0/fixed_step)* sd);
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_sa_H1[i], eps_calc));
       }
    }
    SECTION("Oxygen Small Angle rate"){

      std::vector<double> ref_sa_O16{0.00579506, 0.00763454, 0.011343, 0.0150841, 0.0283483};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        double density = 1.0;
        auto tmp = openmc::proton_sde::moliere_scattering_precomp(O16, E);
        //Single nuclide, so no need to sum anything
        auto sd = openmc::proton_sde::moliere_transform(E, tmp.first, tmp.second, density);
        auto rate =  std::sqrt((1.0/fixed_step)* sd);
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_sa_O16[i], eps_calc));
       }
    }
    SECTION("Carbon Small Angle rate"){

      std::vector<double> ref_sa_C12{0.00741039, 0.00975985, 0.0144957, 0.0192732, 0.036217};
      for(int i = 0; i < 5; i++){
        auto E = energies[i]*1e6;
        double density = 2.0;
        auto tmp = openmc::proton_sde::moliere_scattering_precomp(C12, E);
        //Single nuclide, so no need to sum anything
        auto sd = openmc::proton_sde::moliere_transform(E, tmp.first, tmp.second, density);
        auto rate =  std::sqrt((1.0/fixed_step)* sd);
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_sa_C12[i], eps_calc));
       }
    }

}
 