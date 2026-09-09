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
inline double mock_random_value(){
    double sample = uniform_dist(proton_rng) * 0.1/0.000668456;
    return sample;
}

inline double random_angle(){
  return angle_dist(proton_rng);
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
  const double arbitrary_scale = 1e17;
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
    const double arbitrary_factor = 1e-50;
    return arbitrary_factor * nuclide.proton_ne_rate.evaluate(e/1e6) * nuclide.Z_ / nuclide.A_;
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
    const double arbitrary_factor = 5e1;
    return arbitrary_factor * nuclide.proton_el_rate.evaluate(e/1e6) * nuclide.Z_ / nuclide.A_;
    //return 1e4;
  }

  /* Fortran re-design
    !> \brief Computes the standard deviation of Moliere scattering
    !> Based on p. 7 and 8 of [1]
    !> References:
    !>  [1] https://doi.org/10.1088/1361-6560/ae5586
    !> \param material The material through which the proton travels
    !> \param energy The energy of the proton
    !> \param time_step The time step for the simulation
    !> \return The standard deviation of Moliere scattering
    PURE FUNCTION moliere_scattering_sd(material, energy, time_step) RESULT(sd)
      TYPE(pt_material), INTENT(IN) :: material
      REAL(KIND=REAL64), INTENT(IN) :: energy, time_step
      REAL(KIND=REAL64), PARAMETER :: fixed_time_step = 0.05
      REAL(KIND=REAL64) :: chi_a_sq, chi_c_sq, pv_sq, beta_sq, sd, omega, temp1, temp2
      INTEGER :: i

      ! beta^2 see underneath eq. (3) on p. 6 of [1]
      beta_sq = (2.0 * mpcsq + energy) * energy /(mpcsq + energy)**2

      ! (p*beta)^2
      pv_sq = (2.0 * mpcsq + energy) * energy / (mpcsq + energy)
      pv_sq = pv_sq**2

      chi_c_sq = 0.0_REAL64
      chi_a_sq = 0.0_REAL64
      ! chi_c_sq is sum of individual contributions from each nuclide, 
      ! chi_a_sq is a weighted average on a log scale
      DO i = 1, material%no_nucs
        ! Z_i(Z_i+1)/A_i
        temp1 = material%nucs(i)%massFraction * material%nucs(i)%Z * (material%nucs(i)%Z + 1.0_REAL64) / material%nucs(i)%A
        chi_c_sq = chi_c_sq + temp1
        ! (chi_alpha,i)^2, note pv_sq = (p * beta)^2
        temp2 = 2.007E-5_REAL64 * REAL(material%nucs(i)%Z, KIND=REAL64)**(2.0/3.0) * &
          (1.0_REAL64 + 3.34_REAL64 * (material%nucs(i)%Z * alpha)**2 / beta_sq) * beta_sq / pv_sq
        chi_a_sq = chi_a_sq + temp1 * LOG(temp2)
      END DO
      ! normalise and eliminate the log, 
      ! note denominator in log(chi_a_sq) same as chi_c_sq before multiplying by nucleide independent parameters
      chi_a_sq = EXP(chi_a_sq / chi_c_sq)
      ! multiply chi_c_sq by time_step and and parameters independent of the individual nucleides
      chi_c_sq = chi_c_sq * 0.157_REAL64 * fixed_time_step * material%density / pv_sq
      omega = chi_c_sq / (chi_a_sq * 2.0_REAL64 * (1.0_REAL64 - 0.98_REAL64)) ! 0.98 - truncation parameter, see p. 8 [1]
      ! standard deviation
      sd = SQRT(time_step/fixed_time_step * chi_c_sq * ((1.0_REAL64 + omega) * LOG(1.0_REAL64 + omega) / omega - 1.0_REAL64) / (1.0_REAL64 + 0.98_REAL64**2))
    END FUNCTION
    */
  inline std::pair<double, double> moliere_scattering_precomp(int i_nuclide, double e){
    auto energy = e/1e6; // Converting to MeV
    const double alpha = 1.0/137.0;
    const double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    // beta^2 see underneath eq. (3) on p. 6 of [1]
    auto beta_sq = (2.0 * mpcsq + energy) * energy / std::pow(mpcsq + energy, 2);

    // (p*beta)^2
    auto pv_sq = (2.0 * mpcsq + energy) * energy / (mpcsq + energy);
    pv_sq = pv_sq * pv_sq;

    const Nuclide_t& nuclide = *data::nuclides.at(i_nuclide);
    // chi_c_sq is sum of individual contributions from each nuclide, 
    // chi_a_sq is a weighted average on a log scale
    // HERE we only calculate the per-nuclide part
    //! Z_i(Z_i+1)/A_i
    auto temp1 = nuclide.Z_ * (nuclide.Z_ + 1.0)/nuclide.A_;
        //chi_c_sq = chi_c_sq + temp1
        //! (chi_alpha,i)^2, note pv_sq = (p * beta)^2
    auto temp2 = 2.007E-5 * std::pow(nuclide.Z_, 2.0/3.0) * (1.0 + 3.34 * std::pow(nuclide.Z_ * alpha/beta_sq, 2)) * beta_sq / pv_sq;
    return {temp1, temp1 * log(temp2)};
  }
  inline std::pair<double, double> moliere_transform(double energy, double sum_c, double sum_a, double density){
    const double fixed_step = 0.05;
    auto chi_a_sq = exp(sum_a/sum_c);
    energy = energy / 1e6;
    const double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    auto pv_sq = (2.0 * mpcsq + energy) * energy / (mpcsq + energy);
    pv_sq = pv_sq * pv_sq;
    auto chi_c_sq = sum_c + 0.157 * fixed_step * density / pv_sq;
    auto omega = chi_c_sq / (chi_a_sq * 2.0 * (1.0 - 0.98)); // 0.98
    return {chi_c_sq, omega};
  }
    /*
    ! normalise and eliminate the log, 
      ! note denominator in log(chi_a_sq) same as chi_c_sq before multiplying by nucleide independent parameters
      chi_a_sq = EXP(chi_a_sq / chi_c_sq)
      ! multiply chi_c_sq by time_step and and parameters independent of the individual nucleides
      chi_c_sq = chi_c_sq * 0.157_REAL64 * fixed_time_step * material%density / pv_sq
      omega = chi_c_sq / (chi_a_sq * 2.0_REAL64 * (1.0_REAL64 - 0.98_REAL64)) ! 0.98 - truncation parameter, see p. 8 [1]
      ! standard deviation
      sd = SQRT(time_step/fixed_time_step * chi_c_sq * ((1.0_REAL64 + omega) * LOG(1.0_REAL64 + omega) / omega - 1.0_REAL64) / (1.0_REAL64 + 0.98_REAL64**2))
      */


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

  inline std::pair<double, double> spherical_bm(double distance, double energy, std::vector<double> direction_in, std::pair<double,double> moliere_transformed_precomp){
    /* Re-translated from the Fortran decisions
    !> \brief Simulation of the spherical Brownian motion process
    !> Based on Algorithm 1 in [2] 
    !> References:
    !>  [1] https://doi.org/10.1088/1361-6560/ae5586
    !>  [2] https://doi.org/10.1016/j.spl.2020.108836
    !> \param dt The time step for the simulation
    !> \param energy The energy of the proton
    !> \param material The material through which the proton travels
    !> \param direction_in The current direction of the proton in spherical coordinates
    !> \param state The state of the random number generator
    !> \param b_state The state of the Box-Muller random number generator
    !> \return The new direction of the proton in spherical coordinates
    FUNCTION spherical_bm(dt, energy, material, direction_in) RESULT(direction_out)
      REAL(KIND=REAL64), INTENT(IN) :: dt, energy, direction_in(3)
      TYPE(pt_material), INTENT(IN) :: material
      TYPE(BoxMullerRNGState) :: b_state
      REAL(KIND=REAL64) :: direction_out(2),z(3), u(3), w(3), y, denom, theta
      REAL,EXTERNAL :: random

      b_state%has_cache = .false.
      ! Convert to Cartesian coordinates
      !z(1) = sin(direction_in(1)) * cos(direction_in(2))
      !z(2) = sin(direction_in(1)) * sin(direction_in(2))
      !z(3) = cos(direction_in(1))
      z=direction_in
      y = wright_fisher_diffusion(moliere_scattering_sd(material, energy, dt)**2, b_state)
      theta = 2 * PI * random(1)
      
      ! Set up defaults for when z is near (0, 0, 1)
      u = [1.0_REAL64/SQRT(2.0_REAL64), &
          1.0_REAL64/SQRT(2.0_REAL64), &
          0.0_REAL64]

      ! u = (e_3 -z)/|e_3 - z| (line 3, Algorithm 1 [2])
      denom = SQRT(z(1)**2 + z(2)**2 + (z(3) - 1.0_REAL64)**2)
      u = -z / denom
      u(3) = u(3) + 1.0_REAL64/denom

      !Evaluate expression to the right of O(z) in line 4, Algorithm 1 [2]
      z(1) = 2 * SQRT(y * (1.0_REAL64 - y)) * cos(theta)
      z(2) = 2 * SQRT(y * (1.0_REAL64 - y)) * sin(theta)
      z(3) = 1.0_REAL64 - 2 * y

      ! Evaluate O(z)z with O(z)=I-2uu^T (line 3/4 Algorithm 1 [2])
      w = z - 2 * u * (DOT_PRODUCT(u, z))
      
      ! New direction in spherical coordinates
      direction_out(1) = ACOS(w(3))
      direction_out(2) = ATAN2(w(2), w(1))

      direction_out(1) = dot_product(direction_in,w)/sqrt(dot_product(w,w)) ! might not need the denominator if w unity vector

    END FUNCTION
    */
    std::vector<double> z, u, w;
    u.resize(3);
    w.resize(3);

    z=direction_in;

    auto chi_c_sq = moliere_transformed_precomp.first;
    auto omega = moliere_transformed_precomp.second;
    auto moliere_sd_sq = (distance/0.05 * chi_c_sq * ((1.0 + omega) * log(1.0 + omega) / omega - 1.0) / (1.0 + std::pow(0.98, 2)));
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

      direction_out_1 = dot_product(direction_in,w)/sqrt(dot_product(w,w));
      //! might not need the denominator if w unity vector
      return{direction_out_1, direction_out_2};
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
    nuclide.proton_el_xsec = CS_2d(filename, 0.04);
  }

  inline double rutherford_elastic_scatter(int i_nuclide, double e){
    const Nuclide_t& nuclide = *data::nuclides.at(i_nuclide);
    //std::cout<<"sampling rutherford "<<std::endl;
    auto alpha = nuclide.proton_el_xsec.sample(e/1e6, next_rand());
    //std::cout<<"Done sampling "<<std::endl;
    return cos(alpha); //TODO URGENT cm to lab??
    //return random_angle(); //TODO actual
    //Need the 2D cross section here.
  }

  inline std::pair<double, double> non_elastic_scatter(double initial_energy){
    return {initial_energy - next_rand()*1e3, random_angle()}; // TODO actual
  }

};
#endif