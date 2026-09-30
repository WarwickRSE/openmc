#include <cmath>
#include <vector>
#include <map>

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

// Producing reproducible random sequences for an exact cross-check
// We define one sequence for each test case, and select between them
// using the 'seed' as an identifier.
// This lets us test with more than one sequence
static inline std::map<uint64_t, std::vector<double> > sequences;
static inline std::map<uint64_t, size_t> sequence_index;

inline void register_sequence(uint64_t seed, std::vector<double> & seq){
  if(sequences.count(seed) == 0){
    //register
    sequences[seed] = seq;
    sequence_index[seed] = 0;
  }else{
    throw std::runtime_error("Sequence already registered for seed");
  }
}

inline double prn(uint64_t * seed){
  if(sequences.count(*seed)){
    // Starting the yield
    size_t ind = sequence_index[*seed];
    if(ind < sequences[*seed].size()){
      return sequences[*seed][sequence_index[*seed] ++];
    }else{
      throw std::runtime_error("Random sequence exhausted");
    }
  }else{
    throw std::runtime_error("Random sequence for seed not found");
  }
}
#include "openmc/protons.h"
//------------------------ Basic rate functions ------------------------------------------------

// Threshold for equality, keeping in mind that we might be using
// slightly different calculation breakdowns so do NOT expect full
// FP equality
constexpr double eps_calc = 1e-5;
constexpr double eps_weak = 1e-4;

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
      for(int i = 0; i < energies.size(); i++){
        auto E = energies[i]*1e6;
        auto I = 75.0; //Correct for e.g. H in Water
        auto loss = openmc::proton_sde::proton_bethe_bloch(H1, E, I)/1e6;
        REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_loss_H1[i], eps_calc));
       }
    }
    SECTION("Oxygen Bethe Bloch"){
      // Created by test code
      std::vector<double> ref_loss_O16{4.04867, 4.90762, 6.57085, 8.17145, 13.4217};
      for(int i = 0; i < energies.size(); i++){
        auto E = energies[i]*1e6;
        auto I = 75.0; //Correct for e.g. H in Water
        // Divide by A_ as we canclled this factor
        auto loss = openmc::proton_sde::proton_bethe_bloch(O16, E, I)/1e6/O16.A_;
        REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_loss_O16[i], eps_calc));
      }
    }
    SECTION("Carbon Bethe Bloch"){
      std::vector<double> ref_loss_C12{3.98478, 4.82798, 6.45965, 8.02871, 13.1688};
      for(int i = 0; i < energies.size(); i++){
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
    /* Cross-match to proton_energy_straggle in physics.cpp:
      return std::sqrt(p.macro_xs().energy_straggling * proton_sde::energy_straggling_update_sq(p.E()) * distance) * normal_variate(0.0, 1.0, p.current_seed()) * proton_sde::MeVToeV;
      where p.macro_xs().energy_straggling = micro.energy_straggling;
      and micro.energy_straggling = proton_sde::energy_straggling_sd(nuclide);
    */
    SECTION("Hydrogen Energy Straggling"){
        // Created by test code
      std::vector<double> ref_strag_H1{0.438768, 0.427397, 0.416248, 0.410761, 0.403185};
      for(int i = 0; i < energies.size(); i++){
        auto E = energies[i]*1e6;
        double distance = 1.0;
        auto loss = std::sqrt(openmc::proton_sde::energy_straggling_sd(H1) * openmc::proton_sde::energy_straggling_update_sq(E) * distance);
        REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_strag_H1[i], eps_weak));
       }
    }
    SECTION("Oxygen Energy Straggling"){
        // Created by test code
      std::vector<double> ref_strag_O16{0.311504, 0.303432, 0.295516, 0.291621, 0.286242};
      for(int i = 0; i < energies.size(); i++){
        auto E = energies[i]*1e6;
        double distance = 1.0;
        auto loss = std::sqrt(openmc::proton_sde::energy_straggling_sd(O16) * openmc::proton_sde::energy_straggling_update_sq(E) * distance);
        REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_strag_O16[i], eps_weak));
       }
    }
    SECTION("Carbon Energy Straggling"){
        // Created by test code
      std::vector<double> ref_strag_C12{0.311352, 0.303283, 0.295371, 0.291478, 0.286102};
      for(int i = 0; i < energies.size(); i++){
        auto E = energies[i]*1e6;
        double distance = 1.0;
        auto loss = std::sqrt(openmc::proton_sde::energy_straggling_sd(C12) * openmc::proton_sde::energy_straggling_update_sq(E) * distance);
        REQUIRE_THAT(loss, Catch::Matchers::WithinRel(ref_strag_C12[i], eps_weak));
       }
    }

}

TEST_CASE("Proton Cross Sections from File"){
    //Checking partial contributions from selected Nuclides

    // PRECISE matching here depends on correcting Avogadro's number to 6.022 in the test code to generate comparisons

    const std::string data_path = std::getenv("OPENMC_PROTON_DATA");

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
      std::vector<double> ref_el_H1{0.0400357, 0.0400357, 0.0462034, 0.0589672, 0.11799};
      for(int i = 0; i < energies.size(); i++){
        auto E = energies[i]*1e6;
        double atom_density = 1.0 * openmc::N_AVOGADRO / H1.A_;
        auto rate = openmc::proton_sde::rutherford_elastic_rate(H1, E) * atom_density;
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_el_H1[i], eps_calc));
       }
    }
    SECTION("Oxygen Elastic Rate"){
        // Created by test code
      std::vector<double> ref_el_O16{0.00704476, 0.00704476, 0.0153057, 0.0259321, 0.0765291};
      for(int i = 0; i < energies.size(); i++){
        auto E = energies[i]*1e6;
        double atom_density = 1.0 * openmc::N_AVOGADRO / O16.A_;
        auto rate = openmc::proton_sde::rutherford_elastic_rate(O16, E) * atom_density;
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_el_O16[i], eps_calc));
       }
    }
    SECTION("Carbon Elastic Rate"){
        // Created by test code
      std::vector<double> ref_el_C12{0.013668, 0.013668, 0.0292795, 0.0486201, 0.131744};
      for(int i = 0; i < energies.size(); i++){
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
      std::vector<double> ref_ne_O16{0.011104, 0.011104, 0.0111793, 0.012158, 0.0169007};
      for(int i = 0; i < energies.size(); i++){
        auto E = energies[i]*1e6;
        double atom_density = 1.0 * openmc::N_AVOGADRO / O16.A_;
        auto rate = openmc::proton_sde::non_elastic_rate(O16, E) * atom_density;
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_ne_O16[i], eps_calc));
       }
    }
    SECTION("Carbon Non-Elastic Rate"){
        // Created by test code
      std::vector<double> ref_ne_C12{0.0222966, 0.0222966, 0.0227629, 0.0269244, 0.0364004};
      for(int i = 0; i < energies.size(); i++){
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

    SECTION("Hydrogen Small Angle rate"){

      std::vector<double> ref_sa_H1{0.00394011, 0.00519039, 0.00771194, 0.0102573, 0.019297};
      for(int i = 0; i < energies.size(); i++){
        auto E = energies[i]*1e6;
        double density = 1.0;
        auto tmp = openmc::proton_sde::moliere_scattering_precomp(H1, E);
        auto tmp2 = std::make_tuple(tmp.first, tmp.second, density);
        //Single nuclide, so no need to sum anything
        auto sd = openmc::proton_sde::moliere_transform(E, tmp2, 1.0);
        auto rate =  std::sqrt(sd);
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_sa_H1[i], eps_calc));
       }
    }
    SECTION("Oxygen Small Angle rate"){

      std::vector<double> ref_sa_O16{0.00579506, 0.00763454, 0.011343, 0.0150841, 0.0283483};
      for(int i = 0; i < energies.size(); i++){
        auto E = energies[i]*1e6;
        double density = 1.0;
        auto tmp = openmc::proton_sde::moliere_scattering_precomp(O16, E);
        auto tmp2 = std::make_tuple(tmp.first, tmp.second, density);
        //Single nuclide, so no need to sum anything
        auto sd = openmc::proton_sde::moliere_transform(E, tmp2, 1.0);
        auto rate =  std::sqrt(sd);
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_sa_O16[i], eps_calc));
       }
    }
    SECTION("Carbon Small Angle rate"){

      std::vector<double> ref_sa_C12{0.00741039, 0.00975985, 0.0144957, 0.0192732, 0.036217};
      for(int i = 0; i < energies.size(); i++){
        auto E = energies[i]*1e6;
        double density = 2.0;
        auto tmp = openmc::proton_sde::moliere_scattering_precomp(C12, E);
        auto tmp2 = std::make_tuple(tmp.first, tmp.second, density);
        //Single nuclide, so no need to sum anything
        auto sd = openmc::proton_sde::moliere_transform(E, tmp2, 1.0);
        auto rate =  std::sqrt(sd);
        REQUIRE_THAT(rate, Catch::Matchers::WithinRel(ref_sa_C12[i], eps_calc));
       }
    }

}

// Can we test the Spherical walk sensibly here?

//Now test the cross section sampling
TEST_CASE("Proton Large Angle Scattering"){

    const std::string data_path = std::getenv("OPENMC_PROTON_DATA");

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

    SECTION("Hydrogen elastic sampling"){
      std::uint64_t seed = 1234;
      std::vector<double> random_seq{0.1, 0.1, 0.5, 0.0, 1.0, 0.0, 1.0, 0.1, 0.67};
      register_sequence(seed, random_seq);

      std::vector<double> energies{100.0, 73.0, 5.3, 1.0, 1.0, 160.0, 160.0, 150.0, 150.0};

      std::vector<double> ref_el_ang_H1{0.958581, 0.962766, 0.37857, 1.11851, 0.434487, 1.10332, 0.04, 0.962499, 0.308092};
      //{1.44996, 1.53925, 1.53162, 1.55071, 0.04, 1.54993, 0.04, 1.37805, 0.540784};
      for(size_t i = 0; i < energies.size(); i++){
        auto alpha = H1.proton_el_xsec.sample(energies[i], prn(&seed));
        REQUIRE_THAT(alpha, Catch::Matchers::WithinRel(ref_el_ang_H1[i], eps_calc));
      }
    }
    SECTION("Oxygen elastic sampling"){
      std::uint64_t seed = 2345;
      std::vector<double> random_seq{0.1, 0.1, 0.5, 0.0, 1.0, 0.0, 1.0, 0.1, 0.67};
      register_sequence(seed, random_seq);

      std::vector<double> energies{100.0, 73.0, 5.3, 1.0, 1.0, 160.0, 160.0, 150.0, 150.0};

      std::vector<double> ref_el_ang_O16{0.335586, 0.376031, 0.163273, 3.14159, 0.601333, 3.14159, 0.04, 0.265345, 0.0606393};
      //{0.328358, 0.374877, 0.0564652, 3.14159, 0.04, 3.14159, 0.04, 0.260273, 0.0592945};
      for(size_t i = 0; i < energies.size(); i++){
        auto alpha = O16.proton_el_xsec.sample(energies[i], prn(&seed));
        REQUIRE_THAT(alpha, Catch::Matchers::WithinRel(ref_el_ang_O16[i], eps_calc));
      }
    }
    SECTION("Carbon elastic sampling"){
      std::uint64_t seed = 3456;
      std::vector<double> random_seq{0.1, 0.1, 0.5, 0.0, 1.0, 0.0, 1.0, 0.1, 0.67};
      register_sequence(seed, random_seq);

      std::vector<double> energies{100.0, 73.0, 5.3, 1.0, 1.0, 160.0, 160.0, 150.0, 150.0};

      std::vector<double> ref_el_ang_C12{0.354868, 0.394815, 0.141732, 3.14159, 0.55947, 3.14159, 0.04, 0.276945, 0.0738311};
      //{0.349714, 0.395738, 0.0563128, 3.14159, 0.04, 3.14159, 0.04, 0.276185, 0.0763119};
      for(size_t i = 0; i < energies.size(); i++){
        auto alpha = C12.proton_el_xsec.sample(energies[i], prn(&seed));
        REQUIRE_THAT(alpha, Catch::Matchers::WithinRel(ref_el_ang_C12[i], eps_calc));
      }
    }
    SECTION("Hydrogen elastic angle"){
      std::uint64_t seed = 12341;
      std::vector<double> random_seq{0.1, 0.1, 0.5, 0.0, 1.0, 0.0, 1.0, 0.1, 0.67};
      register_sequence(seed, random_seq);

      std::vector<double> energies{100.0, 73.0, 5.3, 1.0, 1.0, 160.0, 160.0, 150.0, 150.0};

      std::vector<std::vector<double>> ref_el_ang_entire_H1{{0.958581, 32.408}, {0.962766, 23.5972}, {0.37857, 4.58002}, {1.11851, 0.197276}, {0.434487, 0.824133}, {1.10332, 31.3802}, {0.04, 159.725}, {0.962499, 47.2617}, {0.308092, 135.329}};
      for(size_t i = 0; i < energies.size(); i++){
        // NOTE: this impl. has already selected a Nuclide by this point
        // Therefore we should only need 1 random roll per call
        auto vals = openmc::proton_sde::rutherford_elastic_scatter(H1, energies[i]*openmc::proton_sde::MeVToeV, &seed);
        REQUIRE_THAT(acos(vals.second), Catch::Matchers::WithinRel(ref_el_ang_entire_H1[i][0], eps_calc));
        REQUIRE_THAT(vals.first*openmc::proton_sde::eVToMeV, Catch::Matchers::WithinRel(ref_el_ang_entire_H1[i][1], eps_calc));
      }
    }
    SECTION("Oxygen elastic angle"){
      std::uint64_t seed = 23451;
      std::vector<double> random_seq{0.1, 0.1, 0.5, 0.0, 1.0, 0.0, 1.0, 0.1, 0.67};
      register_sequence(seed, random_seq);

      std::vector<double> energies{100.0, 73.0, 5.3, 1.0, 1.0, 160.0, 160.0, 150.0, 150.0};

      std::vector<std::vector<double>> ref_el_ang_entire_O16{{0.335586, 99.2683}, {0.376031, 72.3406}, {0.163273, 5.29117}, {3.14159, 0.778443}, {0.601333, 0.978292}, {3.14159, 122.257}, {0.04, 159.983}, {0.265345, 149.293}, {0.0606393, 149.963}};
      for(size_t i = 0; i < energies.size(); i++){
        // NOTE: this impl. has already selected a Nuclide by this point
        // Therefore we should only need 1 random roll per call
        auto vals = openmc::proton_sde::rutherford_elastic_scatter(O16, energies[i]*openmc::proton_sde::MeVToeV, &seed);
        REQUIRE_THAT(acos(vals.second), Catch::Matchers::WithinRel(ref_el_ang_entire_O16[i][0], eps_calc));
        REQUIRE_THAT(vals.first*openmc::proton_sde::eVToMeV, Catch::Matchers::WithinRel(ref_el_ang_entire_O16[i][1], eps_calc));
      }
    }
    SECTION("Carbon elastic angle"){
      std::uint64_t seed = 34561;
      std::vector<double> random_seq{0.1, 0.1, 0.5, 0.0, 1.0, 0.0, 1.0, 0.1, 0.67};
      register_sequence(seed, random_seq);

      std::vector<double> energies{100.0, 73.0, 5.3, 1.0, 1.0, 160.0, 160.0, 150.0, 150.0};

      std::vector<std::vector<double>> ref_el_ang_entire_C12{{0.354868, 98.9132}, {0.394815, 72.0349}, {0.141732, 5.29113}, {3.14159, 0.716088}, {0.55947, 0.974907}, {3.14159, 111.884}, {0.04, 159.977}, {0.276945, 148.976}, {0.0738311, 149.927}};
      for(size_t i = 0; i < energies.size(); i++){
        auto vals = openmc::proton_sde::rutherford_elastic_scatter(C12, energies[i]*openmc::proton_sde::MeVToeV, &seed);
        REQUIRE_THAT(acos(vals.second), Catch::Matchers::WithinRel(ref_el_ang_entire_C12[i][0], eps_calc));
        REQUIRE_THAT(vals.first*openmc::proton_sde::eVToMeV, Catch::Matchers::WithinRel(ref_el_ang_entire_C12[i][1], eps_calc));
      }
    }

    SECTION("Oxygen non elastic sampling"){
      std::uint64_t seed = 5678;
      std::vector<double> random_seq{0.1, 0.1, 0.5, 0.0, 1.0, 0.0, 1.0, 0.1, 0.67};
      register_sequence(seed, random_seq);

      std::vector<double> energies{100.0, 73.0, 5.3, 1.0, 1.0, 160.0, 160.0, 150.0, 150.0};

      std::vector<std::vector<double>> ref_ne_ang_O16{{0.155904, 1.77733}, {0.214461, 1.62208}, {0, 0}, {0, 0}, {0, 0.1}, {0.04473, 0}, {1, 126.528}, {0.106837, 1.68676}, {0.739901, 27.8417}};
      for(size_t i = 0; i < energies.size(); i++){
        double out_r, out_e;
        O16.proton_ne_xsec.sample(energies[i], out_r, out_e, prn(&seed));
        REQUIRE_THAT(out_r, Catch::Matchers::WithinRel(ref_ne_ang_O16[i][0], eps_calc));
        REQUIRE_THAT(out_e, Catch::Matchers::WithinRel(ref_ne_ang_O16[i][1], eps_calc));
      }
    }
    SECTION("Carbon non elastic sampling"){
      std::uint64_t seed = 6789;
      std::vector<double> random_seq{0.1, 0.1, 0.5, 0.0, 1.0, 0.0, 1.0, 0.1, 0.67};
      register_sequence(seed, random_seq);

      std::vector<double> energies{100.0, 73.0, 5.3, 1.0, 1.0, 160.0, 160.0, 150.0, 150.0};
      
      std::vector<std::vector<double>> ref_ne_ang_C12{{0.171597, 1.40419}, {0.246883, 1.41405}, {0, 0.207567}, {0, 0}, {0, 0.1}, {0.04816, 0}, {1, 123.157}, {0.115362, 1.20864}, {0.925111, 38.1043}};
      for(size_t i = 0; i < energies.size(); i++){
        double out_r, out_e;
        C12.proton_ne_xsec.sample(energies[i], out_r, out_e, prn(&seed));
        REQUIRE_THAT(out_r, Catch::Matchers::WithinRel(ref_ne_ang_C12[i][0], eps_calc));
        REQUIRE_THAT(out_e, Catch::Matchers::WithinRel(ref_ne_ang_C12[i][1], eps_calc));
     }
    }
    SECTION("Oxygen non elastic angle and energy"){
      std::uint64_t seed = 4321;
      std::vector<double> random_seq{0.1, 0.1, 0.1, 0.1, 0.5, 0.5, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.1, 0.1, 0.67, 0.67};
      register_sequence(seed, random_seq);

      std::vector<double> energies{100.0, 73.0, 5.3, 1.0, 1.0, 160.0, 160.0, 150.0, 150.0};

      std::vector<std::vector<double>> ref_ne_ang_entire_O16{{-0.505378, 0.875054}, {-0.545804, 0.860497}, {1, 0}, {1, 0}, {1, 0.140666}, {1, 0}, {1, 143.821}, {-0.375194, 0.711607}, {0.744059, 33.5311}};
      for(size_t i = 0; i < energies.size(); i++){
        auto vals = openmc::proton_sde::non_elastic_scatter(O16, energies[i]*openmc::proton_sde::MeVToeV, &seed);
        //NOTE: items reversed wrt test code
        REQUIRE_THAT(vals.second, Catch::Matchers::WithinRel(ref_ne_ang_entire_O16[i][0], eps_calc));
        REQUIRE_THAT(vals.first*openmc::proton_sde::eVToMeV, Catch::Matchers::WithinRel(ref_ne_ang_entire_O16[i][1], eps_calc));
      }
    }
    SECTION("Carbon non elastic angle and energy"){
      std::uint64_t seed = 5432;
      std::vector<double> random_seq{0.1, 0.1, 0.1, 0.1, 0.5, 0.5, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.1, 0.1, 0.67, 0.67};
      register_sequence(seed, random_seq);

      std::vector<double> energies{100.0, 73.0, 5.3, 1.0, 1.0, 160.0, 160.0, 150.0, 150.0};

      std::vector<std::vector<double>> ref_ne_ang_entire_C12{{-0.234352, 0.547033}, {-0.361884, 0.611241}, {0.362028, 0.238875}, {1, 0}, {1, 0.154516}, {1, 0}, {1, 145.68}, {0.0962637, 0.443222}, {0.83307, 48.0949}};
      for(size_t i = 0; i < energies.size(); i++){
        auto vals = openmc::proton_sde::non_elastic_scatter(C12, energies[i]*openmc::proton_sde::MeVToeV, &seed);
        //NOTE: items reversed wrt test code
        REQUIRE_THAT(vals.second, Catch::Matchers::WithinRel(ref_ne_ang_entire_C12[i][0], eps_calc));
        REQUIRE_THAT(vals.first*openmc::proton_sde::eVToMeV, Catch::Matchers::WithinRel(ref_ne_ang_entire_C12[i][1], eps_calc));
      }
    }
 

}