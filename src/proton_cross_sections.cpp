#include <unordered_map>

#include "openmc/proton_cross_sections.h"
#include "openmc/nuclide.h"

namespace openmc{
  //TODO - writen by AI because lazy but better to name files using chemical symbols in the first place??
  std::string letter_to_string(std::string sym){
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

  void read_proton_data(int i_nuclide, std::string name){
    //NOTE: reads for EACH isotope afresh. TODO - fix...
    std::string path = "/media/raid/MathRadData/protons/";
    std::string filename = path;
    auto pos = name.find_first_of("0123456789");
    //TODO - better....
    auto sym = letter_to_string(name.substr(0, pos));
    filename += sym;
    filename += "_ne_rate.txt";
    std::cout<<filename<<std::endl;
    Nuclide& nuclide = *data::nuclides.at(i_nuclide);
    nuclide.proton_ne_rate = CS_1d(filename);
    nuclide.proton_ne_rate.check();
    // TODO remove double read
    filename = path + sym+ "_el_ruth_cross_sec.txt";
     std::cout<<filename<<std::endl;
    nuclide.proton_el_rate = CS_1d(filename);
    nuclide.proton_el_rate.check();
    nuclide.proton_el_xsec = CS_2d(filename, 0.04);

    filename = path + sym+ "_ne_energyangle_cdf.txt";
    nuclide.proton_ne_xsec = CS_3d(filename);
  }

}