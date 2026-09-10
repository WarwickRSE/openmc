#ifndef __proton_physics__
#define __proton_physics__

#include <random>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
//TODO - what RNG to use?? NOT a file local static one pls!!
// At that point also move the distribs into the relavant functions

#include "openmc/nuclide.h"
#include "openmc/proton_cross_sections.h"

namespace openmc{
static inline std::mt19937 proton_rng {std::random_device {}()};
constexpr double MAX_DEFLECTION = 1.0e-2; // radians
constexpr double fixed_step = 0.05; // A fixed step length used in temporary calc
constexpr double MeVToeV = 1e6;
constexpr double eVToMeV = 1e-6;
constexpr double alpha_finestruc = 1.0/137.0;
constexpr double mecsq = 0.511;   // mass of electron * speed of light squared, MeV
constexpr double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
constexpr double log_hbar = -21 * log(10) + log(4.136) - log(2 * M_PI); // MeV * s
constexpr double log_c = log(29979245800);                              // cm / s
constexpr double log_avogadro = log(6) + 23 * log(10);
 

static inline std::uniform_real_distribution<double> uniform_dist {0.0, 1.0};
static inline std::uniform_real_distribution<double> angle_dist {std::cos(MAX_DEFLECTION), 1.0};
static inline std::normal_distribution<double> e_strag(0.0, 1.0);
static inline std::normal_distribution<double> generic_gauss(0.0, 1.0);

// IMPORTANT - THIS IS A WIP. A lot of this file is dumb static global state in order to test the MODEL needs before integrating to the codebase proper
//TODO - move some of this into the Nuclide, Material or Particle classes?

inline double log_beta_fn(int a_in, int b_in){
  double a = (double)a_in;
  double b = (double)b_in;
  return std::lgamma(a) + std::lgamma(b) - std::lgamma(a+b);
}

inline double sample_beta(int beta){
  auto ran = log_beta_fn(1+beta, 1)+log((1.0+beta) * uniform_dist(proton_rng));
  ran = ran*(1.0/(1.0 + beta));
  ran = 1.0 - exp(ran);
  return ran;
}

inline double log_factorial(int n) {
  return std::lgamma(static_cast<double>(n) + 1.0);
}

inline double log_pochhammer(double a, double x){
    return std::lgamma(a + x) - std::lgamma(a);
}

inline double next_rand(){
  return uniform_dist(proton_rng);
}

inline double random_exp(double lambda){
  std::exponential_distribution<double> generic_exp(lambda);
  return generic_exp(proton_rng);
}

  inline double dot_product(std::vector<double> a, std::vector<double> b){
    //Dot product IFF length of a and b is 3
    double val = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    return val;
  }

  /** @brief Calculate BetaSq factor
   * 
   * Used by many of the calculations, See Eq (3), p 6 of [1].
   * @param E Energy in MeV
   */
  inline constexpr double betaSq(double E){return (2.0 * mpcsq + E) * E / pow(mpcsq + E, 2);}
  /**
   * @brief Calculate pv_sq
   * 
   * This is (p*beta)**2
   * @param E Energy in MeV
   */
  inline constexpr double pvSq(double E){auto tmp = (2.0 * mpcsq + E) * E / (mpcsq + E); return tmp * tmp;}

/** @brief Inelastic energy loss
 * 
 * Computes the inelastic energy loss per cm using bethe-bloch formula
 * This is for a single nuclide and is per density. Multiply by density to get an energy loss in eV/cm
 * 
 * @param i_nuclide The index for this nuclide in the global table
 * @param E Initial Energy of the proton in eV
 * @param I Mean activation energy for current material in MeV // TODO pass as eV
 */
inline double proton_bethe_bloch(int i_nuclide, double E, double I){
    const Nuclide& nuclide = *data::nuclides.at(i_nuclide);
    // Bethe-bloch contrib for single atom of THIS Nuclide only, summed later
    E = E * eVToMeV; // Inside here, expecting MeV
    double betasq = betaSq(E);
    return MeVToeV * 0.3072 * nuclide.Z_ *
             (log(2 * mecsq * betasq / (I * (1 - betasq))) - betasq) /
             (betasq);
}

inline double random_straggle(){
  return e_strag(proton_rng);
}

//Random bump - Energy loss or gain due to straggling. Depends on e and material
// and introduced a random Gaussian
inline double energy_straggling_update_sq(double e){
    e = e*eVToMeV;
    double betasq = betaSq(e);
    return 4 * PI * (1 - betasq / 2) / (1 - betasq) *
        exp(2 * (log(alpha_finestruc) + log_hbar + log_c));
}
  //Duplicated from SDE code
  inline double energy_straggling_sd(int i_nuclide) {

    //The micro part is just the sum of x * z / a; and can be cached, electrons per average molecule in this material
    //The REST is based on the energy
    const Nuclide& nuclide = *data::nuclides.at(i_nuclide);
    return nuclide.Z_ / nuclide.A_;
  }

  inline double non_elastic_rate(int i_nuclide, double e){
    const Nuclide& nuclide = *data::nuclides.at(i_nuclide);
    return nuclide.proton_ne_rate.evaluate(e*eVToMeV);
  }

   /** @brief Retrieve scattering rate for Rutherford and elastic scattering
   * 
   * @param e energy of particle being scattered in eV
   * @return  scattering rate sigma_e 
   */
  inline double rutherford_elastic_rate(int i_nuclide, double e){
    const Nuclide& nuclide = *data::nuclides.at(i_nuclide);
    return nuclide.proton_el_rate.evaluate(e*eVToMeV);
  }

  inline std::pair<double, double> moliere_scattering_precomp(int i_nuclide, double e){
    auto energy = e*eVToMeV;
    auto beta_sq = betaSq(e);
    auto pv_sq = pvSq(e);

    const Nuclide& nuclide = *data::nuclides.at(i_nuclide);
    // chi_c_sq is sum of individual contributions from each nuclide, 
    // chi_a_sq is a weighted average on a log scale
    // HERE we only calculate the per-nuclide part
    //! Z_i(Z_i+1)/A_i
    auto temp1 = nuclide.Z_ * (nuclide.Z_ + 1.0)/nuclide.A_;
        //chi_c_sq = chi_c_sq + temp1
        //! (chi_alpha,i)^2, note pv_sq = (p * beta)^2
    auto temp2 = 2.007E-5 * std::pow(nuclide.Z_, 2.0/3.0) * (1.0 + 3.34 * std::pow(nuclide.Z_ * alpha_finestruc, 2)/beta_sq) * beta_sq / pv_sq;
    return {temp1, temp1 * log(temp2)};
  }
  inline double moliere_transform(double energy, double sum_c, double sum_a, double density){
    auto chi_a_sq = exp(sum_a/sum_c);
    energy = energy * eVToMeV;
    auto pv_sq = pvSq(energy);
    auto chi_c_sq = sum_c * 0.157 * fixed_step * density / pv_sq;
    auto omega = chi_c_sq / (chi_a_sq * 2.0 * (1.0 - 0.98)); // 0.98
    return chi_c_sq * ((1.0 + omega) * log(1.0 + omega) / omega - 1.0) / (1.0 + std::pow(0.98, 2));
  }

  //From the test code

  inline double log_a(const int k, const int m, const double theta){
    double ret = log(theta + 2 * k - 1) + log_pochhammer(theta + m, k - 1) -
                 log_factorial(m) - log_factorial(k - m);
    return ret;
  }

  inline double b(const int k, const int m, const double t, const double theta) {
    double ret = 1;
    if (k > 0) {
      ret = exp(log_a(k, m, theta) - k * (k + theta - 1) * t / 2);
    }
    return ret;
  }

  inline int c(const int m, const double t, const double theta){
    int i = 0;
    double b_curr = b(m, m, t, theta);
    double b_next = b(m + 1, m, t, theta);
    while (b_next >= b_curr) {
      i++;
      b_curr = b_next;
      b_next = b(i + m + 1, m, t, theta);
    }
    return i;
  }

  inline int number_of_blocks(double t) {
    int m = 0;
    double theta = 1;
    if (t < 0.07) {
      double mu = 2.0 / t;
      double sigma = sqrt(2.0 / (3.0 * t));
      m = round(mu + sigma * generic_gauss(proton_rng));
      //m = round(mu + sigma * gsl_ran_gaussian_ziggurat(gen, 1));
    } else {
      std::vector<int> k(1, 0);
      bool proceed = true;
      double u = 0.0;
      while(u == 0.0){
        u = uniform_dist(proton_rng); // TODO URGENT skip 0 somehow less stupid
      }
      double smin = 0, smax = 0, increment = 0;
      while (proceed) {
        k[m] = ceil(c(m, t, theta) / 2.0);
        for (int i = 0; i < k[m]; i++) {
          increment = b(m + 2 * i, m, t, theta) - b(m + 2 * i + 1, m, t, theta);
          smin += increment;
          smax += increment;
        }
        increment = b(m + 2 * k[m], m, t, theta);
        smin += increment - b(m + 2 * k[m] + 1, m, t, theta);
        smax += increment;
        while (smin < u && u < smax) {
          for (int i = 0; i <= m; i++) {
            k[i]++;
            increment = b(i + 2 * k[i], i, t, theta);
            smax = smin + increment;
            smin += increment - b(i + 2 * k[i] + 1, i, t, theta);
          }
        }
        if (smin > u) {
          proceed = false;
        } else {
          k.push_back(0);
          m++;
        }
      }
    }
    return m;
  }

  inline double wright_fisher_diffusion(double r){
    double y;
    if (r > 1e-9) {
      int m = number_of_blocks(r);
      y = sample_beta(1 + m);
    } else {
      y = r / 2;
      std::normal_distribution<double> distribution(0.0, sqrt(r * y * (1 - y)));
      auto sample = distribution(proton_rng);
      y = std::abs(sample);
    }
    return y;
  }

  inline std::pair<double, double> spherical_bm(double distance, double energy, std::vector<double> direction_in, double moliere_transformed_precomp){
    std::vector<double> z, u, w;
    u.resize(3);
    w.resize(3);

    z = direction_in;

    auto moliere_sd_sq = (distance/fixed_step)* moliere_transformed_precomp;
    auto y = wright_fisher_diffusion(moliere_sd_sq);
    auto theta = 2.0 * PI * next_rand();

      // Set up defaults for when z is near (0, 0, 1)
      u = {1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0};

      // u = (e_3 -z)/|e_3 - z| (line 3, Algorithm 1 [2])
      auto denom = sqrt(z[0]*z[0] + z[1]*z[1] + (z[2] - 1.0)*(z[2] - 1.0));
      if (denom > 1e-10) {
        u[0] = -z[0] / denom;
        u[1] = -z[1] / denom;
        u[2] = (1 - z[2]) / denom;
      }

      z[0] = 2 * sqrt(y * (1 - y)) * cos(theta);
      z[1] = 2 * sqrt(y * (1 - y)) * sin(theta);
      z[2] = 1 - 2 * y;
      w[0] = (1 - 2 * u[0] * u[0]) * z[0] - 2 * u[0] * u[1] * z[1] -
           2 * u[0] * u[2] * z[2];
      w[1] = (1 - 2 * u[1] * u[1]) * z[1] - 2 * u[0] * u[1] * z[0] -
           2 * u[1] * u[2] * z[2];
      w[2] = (1 - 2 * u[2] * u[2]) * z[2] - 2 * u[0] * u[2] * z[0] -
           2 * u[1] * u[2] * z[1];

      // New direction in spherical coordinates
      auto direction_out_1 = acos(w[2]);
      auto direction_out_2 = atan2(w[1], w[0]);

      //TODO either actualyl Fake direction in, and skip the extra checks OR pass the real direction and update it
      //direction_out_1 = dot_product(direction_in,w)/sqrt(dot_product(w,w));
      //! might not need the denominator if w unity vector
      return{cos(direction_out_1), direction_out_2};
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

  inline double rutherford_elastic_scatter(int i_nuclide, double e){
    const Nuclide& nuclide = *data::nuclides.at(i_nuclide);
    auto alpha = nuclide.proton_el_xsec.sample(e*eVToMeV, next_rand());
    return cos(alpha); //TODO URGENT cm to lab??
  }

  /** Compute separation energies for incident particles (S_a)
   * S_a as in Section 6.2.3.2 on p. 137 [4]
   * References:
   * [4] https://doi.org/10.2172/1425114
   * @param a Atomic Weight of Nuclide
   * @param z Atomic Number of Nuclide
   * @return S_a
   */
  inline double s(double a, double z) {
    const double a_c = a + 1;
    const double n_c = a - z;
    const double z_c = z + 1;
    const double a_a = a;
    const double n_a = a - z;
    const double z_a = z;
    return 15.68 * (a_c - a_a) -
                 28.07 * (pow(n_c - z_c, 2) / a_c - pow(n_a - z_a, 2) / a_a) -
                 18.56 * (pow(a_c, 2.0 / 3) - pow(a_a, 2.0 / 3)) +
                 33.22 * (pow(n_c - z_c, 2) / pow(a_c, 4.0 / 3) -
                          pow(n_a - z_a, 2) / pow(a_a, 4.0 / 3)) -
                 0.717 * (z_c * z_c / pow(a_c, 1.0 / 3) -
                          z_a * z_a / pow(a_a, 1.0 / 3)) +
                 1.211 * (z_c * z_c / a_c - z_a * z_a / a_a);
  }

   /** Sample from distribution of nonelastic collision
   * 
   * Samples an inelastic collision against a specific nuclide, N
   * References:
   * [1] https://doi.org/10.1088/1361-6560/ae5586
   * [4] https://doi.org/10.2172/1425114
   * @param nuclide The nuclide to collide with, N
   * @param e Energy of incident proton, will be updated
   * @param alpha Polar angle of incident proton, will be updated
   * @param u A uniform random variate in [0,1]
   * @param u2 A uniform random variate in [0,1]
   */
  inline void sample_nonelastic_collision(const openmc::Nuclide & nuclide, double &e, double &alpha, double u, double u2){
    double out_rvalue, out_energy_cm;
    double a = nuclide.A_, z = nuclide.Z_;
    nuclide.proton_ne_xsec.sample(e, out_rvalue, out_energy_cm, u);
    double eps_a = a * e / (a + 1);
    double eps_b = (a + 1) * out_energy_cm / a;
    double e_a = eps_a + s(a, z);
    double e_b = eps_b + s(a, z);
    double x1 = fmin(e_a, 130) * e_b / e_a;
    double x3 = fmin(e_a, 41) * e_b / e_a;
    double aval = 0.04 * x1 + 1.8 * 1e-6 * pow(x1, 3) + 6.7 * 1e-7 * pow(x3, 4);
    double cdfc2 = out_rvalue * cosh(aval) - sinh(aval);
    double cdfc1 = 2 * sinh(aval);
    //double u2 = gsl_rng_uniform(gen);
    double z1 = cdfc1 * u2 + cdfc2;
    double z2 =
        (z1 + sqrt(pow(z1, 2) - pow(out_rvalue, 2) + 1)) / (out_rvalue + 1);
    double out_angle_cm = log(z2) / aval;
    double out_energy_lab =
        out_energy_cm + e / pow(a + 1, 2) +
        2 * sqrt(out_energy_cm * e) * out_angle_cm / (a + 1);
    double out_angle_lab = sqrt(out_energy_cm / out_energy_lab) * out_angle_cm +
                           sqrt(e / out_energy_lab) / (a + 1);
    e = out_energy_lab;
    alpha = out_angle_lab;
    if (out_energy_cm == 0) {
      e = 0;
      alpha = 1; // If outgoing energy is 0, then out_angle_lab should be 1,
                 // rounding errors allow it to be slightly above 1 which is
                 // invalid.
    }
  }

  /** @brief Evaluate an inelastic collision event
   *
   * Evaluates a single inelastic collision with a single Nuclide for a proton at energy e.
   * 
   * @param i_nuclide Index of the nuclide to collide with
   * @param e Energy of the incident proton in eV
   * @return A pair, containing the updated energy in eV and the polar scattering angle
   */
  inline std::pair<double, double> non_elastic_scatter(int i_nuclide, double e){
    const Nuclide& nuclide = *data::nuclides.at(i_nuclide);
    double alpha;
    e = e*eVToMeV; // e passed by value so working with a COPY below
    sample_nonelastic_collision(nuclide, e, alpha, next_rand(), next_rand());
    return {e*MeVToeV, alpha};
  }

};
#endif