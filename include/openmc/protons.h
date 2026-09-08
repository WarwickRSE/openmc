#ifndef __proton_physics__
#define __proton_physics__

#include <random>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
//TODO - what RNG to use?? NOT a file local static one pls!!

#include "openmc/nuclide.h"
using Nuclide_t = openmc::Nuclide;
//#include "fake_nuclide.h"

#include "openmc/proton_cross_sections.h"

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

  inline double non_elastic_rate(int i_nuclide, double e){
    //double log_avogadro = log(6) + 23 * log(10);
    //double log_barns_to_cmsq = -24 * log(10);
    //double ret = 0;
    const Nuclide_t& nuclide = *data::nuclides.at(i_nuclide);
    return nuclide.proton_ne_rate.evaluate(e) * nuclide.Z_ / nuclide.A_;
    //TODO URGENT units
    /*for (unsigned int i = 0; i < at.size(); i++) {
      ret += x[i] * at[i].ne_rate.evaluate(e) / at[i].a;
    }
    double log_molecule_density =
        log(density) + log_avogadro; // molecules / cm^3
    //ret *= exp(log_barns_to_cmsq + log_molecule_density);
    return ret; // rate per cm*/
  }

  inline double rutherford_elastic_rate(int i_nuclide, double e){
    /*double log_avogadro = log(6) + 23 * log(10);
    double log_barns_to_cmsq = -24 * log(10);
    double ret = 0;
    for (unsigned int i = 0; i < at.size(); i++) {
      ret += x[i] * at[i].el_ruth_rate.evaluate(e) / at[i].a;
    }
    double log_molecule_density =
        log(density) + log_avogadro; // molecules / cm^3
    ret *= exp(log_barns_to_cmsq + log_molecule_density);
    return ret; // rate per cm*/
    const Nuclide_t& nuclide = *data::nuclides.at(i_nuclide);
    return nuclide.proton_el_rate.evaluate(e) * nuclide.Z_ / nuclide.A_;
    //return 1e4;
  }

  //TODO - writen by AI because lazy but better to name files using chemical symbols in the first place??
  inline std::string letter_to_string(std::string sym){
    static const std::unordered_map<std::string, std::string> names {
      {"H", "hydrogen"}, {"He", "helium"}, {"Li", "lithium"},
      {"Be", "beryllium"}, {"B", "boron"}, {"C", "carbon"},
      {"N", "nitrogen"}, {"O", "oxygen"}, {"F", "fluorine"},
      {"Ne", "neon"}, {"Na", "sodium"}, {"Mg", "magnesium"},
      {"Al", "aluminum"}, {"Si", "silicon"}, {"P", "phosphorus"},
      {"S", "sulfur"}, {"Cl", "chlorine"}, {"Ar", "argon"},
      {"K", "potassium"}, {"Ca", "calcium"}, {"Sc", "scandium"},
      {"Ti", "titanium"}, {"V", "vanadium"}, {"Cr", "chromium"},
      {"Mn", "manganese"}, {"Fe", "iron"}, {"Co", "cobalt"},
      {"Ni", "nickel"}, {"Cu", "copper"}, {"Zn", "zinc"},
      {"Ga", "gallium"}, {"Ge", "germanium"}, {"As", "arsenic"},
      {"Se", "selenium"}, {"Br", "bromine"}, {"Kr", "krypton"},
      {"Rb", "rubidium"}, {"Sr", "strontium"}, {"Y", "yttrium"},
      {"Zr", "zirconium"}, {"Nb", "niobium"}, {"Mo", "molybdenum"},
      {"Tc", "technetium"}, {"Ru", "ruthenium"}, {"Rh", "rhodium"},
      {"Pd", "palladium"}, {"Ag", "silver"}, {"Cd", "cadmium"},
      {"In", "indium"}, {"Sn", "tin"}, {"Sb", "antimony"},
      {"Te", "tellurium"}, {"I", "iodine"}, {"Xe", "xenon"},
      {"Cs", "cesium"}, {"Ba", "barium"}, {"La", "lanthanum"},
      {"Ce", "cerium"}, {"Pr", "praseodymium"}, {"Nd", "neodymium"},
      {"Pm", "promethium"}, {"Sm", "samarium"}, {"Eu", "europium"},
      {"Gd", "gadolinium"}, {"Tb", "terbium"}, {"Dy", "dysprosium"},
      {"Ho", "holmium"}, {"Er", "erbium"}, {"Tm", "thulium"},
      {"Yb", "ytterbium"}, {"Lu", "lutetium"}, {"Hf", "hafnium"},
      {"Ta", "tantalum"}, {"W", "tungsten"}, {"Re", "rhenium"},
      {"Os", "osmium"}, {"Ir", "iridium"}, {"Pt", "platinum"},
      {"Au", "gold"}, {"Hg", "mercury"}, {"Tl", "thallium"},
      {"Pb", "lead"}, {"Bi", "bismuth"}, {"Po", "polonium"},
      {"At", "astatine"}, {"Rn", "radon"}, {"Fr", "francium"},
      {"Ra", "radium"}, {"Ac", "actinium"}, {"Th", "thorium"},
      {"Pa", "protactinium"}, {"U", "uranium"}, {"Np", "neptunium"},
      {"Pu", "plutonium"}, {"Am", "americium"}, {"Cm", "curium"},
      {"Bk", "berkelium"}, {"Cf", "californium"}, {"Es", "einsteinium"},
      {"Fm", "fermium"}, {"Md", "mendelevium"}, {"No", "nobelium"},
      {"Lr", "lawrencium"}, {"Rf", "rutherfordium"}, {"Db", "dubnium"},
      {"Sg", "seaborgium"}, {"Bh", "bohrium"}, {"Hs", "hassium"},
      {"Mt", "meitnerium"}, {"Ds", "darmstadtium"}, {"Rg", "roentgenium"},
      {"Cn", "copernicium"}, {"Nh", "nihonium"}, {"Fl", "flerovium"},
      {"Mc", "moscovium"}, {"Lv", "livermorium"}, {"Ts", "tennessine"},
      {"Og", "oganesson"}
    };

    auto it = names.find(sym);
    if (it == names.end()) {
      throw std::invalid_argument("Unknown element symbol: " + sym);
    }
    return it->second;
  }

  inline void read_proton_data(int i_nuclide, std::string name){
    //NOTE: reads for EACH isotope afresh. TODO - fix...
    std::string path = "/media/raid/MathRadData/protons/";
    std::string filename = path;
    auto pos = name.find_first_of("0123456789");
    //TODO - better....
    auto sym = letter_to_string(name.substr(0, pos));
    filename += sym;
    filename += "_ne_rate.txt";
    std::cout<<filename<<std::endl;
    Nuclide_t& nuclide = *data::nuclides.at(i_nuclide);
    nuclide.proton_ne_rate = CS_1d(filename);
    nuclide.proton_ne_rate.check();
    // TODO remove double read
    filename = path + sym+ "_el_ruth_cross_sec.txt";
     std::cout<<filename<<std::endl;
    nuclide.proton_el_rate = CS_1d(filename);
    nuclide.proton_el_rate.check();
  }

};
#endif