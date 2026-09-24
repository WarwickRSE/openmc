#ifndef __proton_physics__
#define __proton_physics__

#include <random>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

#include "openmc/math_functions.h"
#include "openmc/random_dist.h"

// Allows mocking out of nuclide
#ifndef __nuclide_included__
#include "openmc/nuclide.h"
#define __nuclide_included__
#endif
#include "openmc/proton_cross_sections.h"

namespace openmc{
  
  namespace proton_sde{
    constexpr double MeVToeV = 1e6;
    constexpr double eVToMeV = 1e-6;
    constexpr double alpha_finestruc = 1.0/137.0;
    constexpr double mecsq = 0.511;   // mass of electron * speed of light squared, MeV
    constexpr double mpcsq = 938.346; // mass of proton * speed of light squared, MeV
    // The following constants are evaluated so we can mark constexpr and use in constexpr contexts
    // but before C++26 the log function is not marked as such
    constexpr double log_hbar = -48.7724349; // -21 * log(10) + log(4.136) - log(2 * PI); // MeV * s
    constexpr double log_c = 24.1237712; //log(29979245800);// cm / s

    //TODO - move some of this into the Nuclide, Material or Particle classes?

    //---------Helper functions-------------------------------------
    // Functions for calculating factors used in many functions. Makes it cleaner to read

    //! Calculate BetaSq factor
    //!  
    //! Beta-sq - related to total energy. Used by many of the calculations, See Eq (3), p 6 of [1].
    //! \param E Energy in MeV
    inline constexpr double betaSq(double E){return (2.0 * mpcsq + E) * E / pow(mpcsq + E, 2);}
    //! Calculate pv_sq, which is (p*beta)**2
    //!  \param E Energy in MeV
    inline constexpr double pvSq(double E){auto tmp = (2.0 * mpcsq + E) * E / (mpcsq + E); return tmp * tmp;}

    //! Inelastic energy loss
    //! 
    //! Computes the inelastic energy loss per cm using bethe-bloch formula
    //!  This is for a single nuclide and is per density. Multiply by density to get an energy loss in eV/cm
    //! \param nuclide The Nuclide
    //! \param E Initial Energy of the proton in eV
    //! \param I Mean activation energy for current material in eV
    inline double proton_bethe_bloch(const Nuclide & nuclide, double E, double I){
      //const Nuclide& nuclide = *data::nuclides.at(i_nuclide);
      // Bethe-bloch contrib for single atom of THIS Nuclide only, summed later
      E = E * eVToMeV; // Inside here, expecting MeV
      I = I * eVToMeV; // Ditto
      double betasq = betaSq(E);
      return MeVToeV * 0.3072 * nuclide.Z_ *
             (log(2 * mecsq * betasq / (I * (1 - betasq))) - betasq) /
             (betasq);
    }

    //! Energy straggling prefactor
    //!
    //! Calculate the prefactor for energy straggling adjustment. Depends on energy of incident particle
    //! \param e Energy of incident particle (eV)
    inline constexpr double energy_straggling_update_sq(double e){
      e = e*eVToMeV;
      double betasq = betaSq(e);
      return 4 * PI * (1 - betasq / 2) / (1 - betasq) *
        exp(2 * (log(alpha_finestruc) + log_hbar + log_c)) *6.022e23;
    }
    //! Per nuclide contribution to energy straggling
    //!
    //! Per-nuclide contribution to energy straggling adjustment
    //! \param nuclide The Nuclide
    inline double energy_straggling_sd(const Nuclide & nuclide) {

    //The micro part is just the sum of z; and can be cached, electrons per average molecule in this material
    //The REST is based on the energy
    // TODO - rho_i/A_i is the mass fraction? Need to nail this down
    return nuclide.Z_/nuclide.A_;
  }

    //! Scattering rate for inelastic scattering
    //!
    //! Per nuclide contribution to inelastic scatterint rate
    //! \param nuclide The Nuclide
    //! \param e energy of particle being scattered in eV
    //! \return  The scattering rate sigma_e
    inline double non_elastic_rate(const Nuclide & nuclide, double e){
      return nuclide.proton_ne_rate.evaluate(e*eVToMeV);
    }
    //! Scattering rate for Rutherford and elastic scattering
    //!
    //! Per nuclide contribution to elastic scattering rate 
    //! \param nuclide The Nuclide
    //! \param e energy of particle being scattered in eV
    //! \return  The scattering rate sigma_e
    inline double rutherford_elastic_rate(const Nuclide & nuclide, double e){
      return nuclide.proton_el_rate.evaluate(e*eVToMeV);
    }

    //! Computes the partial nuclide dependent factors in moliere scattering
    //! 
    //! Effectively computes chi_c and chi_a from pages 7 and 8 of [1] for a single species.
    //! 
    //! NOTE: the return here is not strictly chi_c**2, only the per-nuclide part. We multiply in the other factors in moliere_transform. However this factor is the correct one for computing log(chi_a**2) from page 8 of [1]. 
    //! \param nuclide The Nuclide
    //! \param e The energy of the incident proton in eV
    //! \return A pair, consisting of chi_c**2, and chi_c**2 * log(chi_a**2)
    inline std::pair<double, double> moliere_scattering_precomp(const Nuclide & nuclide, double e){
    auto energy = e*eVToMeV;
    auto beta_sq = betaSq(energy);
    auto pv_sq = pvSq(energy);

    // Nuclide dependent factor part 1
    auto temp1 = nuclide.Z_ * (nuclide.Z_ + 1.0)/nuclide.A_;
    //chi_a**2 for single Nuclide
    auto temp2 = 2.007E-5 * std::pow(nuclide.Z_, 2.0/3.0) * (1.0 + 3.34 * std::pow(nuclide.Z_ * alpha_finestruc, 2)/beta_sq) * beta_sq / pv_sq;
    return {temp1, temp1 * log(temp2)};
  }

    //! Moliere's small angle elastic
    //!
    //! See P 7 and 8 of [1]
    //! chi_c and chi_a and omega are partial factors which are not named.
    //! 0.98 is a truncation parameter F and is fixed here
    //! This computes the answer for a small fixed spatial step. The actual distance must be multiplied in later
    //! \param energy energy of particle being scattered in eV
    //! \param sum_c first part of partial calculation
    //! \param sum_a second part of partial calculation
    //! \param density density of the material
    //! \return squared std-deviation sigma-E for this process for distance step
    inline constexpr double moliere_transform(double energy, std::tuple<double, double, double> moliere_facs, double step){
      auto sum_c = std::get<0>(moliere_facs);
      auto sum_a = std::get<1>(moliere_facs);
      auto density = std::get<2>(moliere_facs);

      auto chi_a_sq = exp(sum_a/sum_c);
      energy = energy * eVToMeV;
      auto pv_sq = pvSq(energy);
      auto chi_c_sq = sum_c * 0.157 * step * density / pv_sq;
      auto omega = chi_c_sq / (chi_a_sq * 2.0 * (1.0 - 0.98)); // 0.98
      return chi_c_sq * ((1.0 + omega) * log(1.0 + omega) / omega - 1.0) / (1.0 + std::pow(0.98, 2));
    }

    //! Factor involved in spherical Brownian Motion 
    inline constexpr double log_a(const int k, const int m, const double theta){
      return log(theta + 2 * k - 1) + log_pochhammer(theta + m, k - 1) -
                 log_factorial(m) - log_factorial(k - m);
    }
    //! Factor involved in spherical Brownian Motion 
    inline constexpr double b(const int k, const int m, const double t, const double theta) {
      return (k > 0) ? ( exp(log_a(k, m, theta) - k * (k + theta - 1) * t / 2) ) : 1;
    }
    //! Factor involved in spherical Brownian Motion 
    inline constexpr int c(const int m, const double t, const double theta){
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
    //! Factor involved in spherical Brownian Motion
    inline int number_of_blocks(double t, std::uint64_t * seed) {
      int m = 0;
      double theta = 1;
      if (t < 0.07) { // TODO citation for this cutoff
        double mu = 2.0 / t;
        double sigma = sqrt(2.0 / (3.0 * t));
        m = round(mu + sigma * normal_variate(0.0, 1.0, seed));
      } else {
        std::vector<int> k(1, 0);
        bool proceed = true;
        double u = 0.0;
        while(u == 0.0){
          // Redrawing if we get exact 0.0 - this is the simplest solution which preserves the random qualities
          u = prn(seed);
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

    //! Evaluate Wright Fisher Diffusion process
    //!
    //!
    inline double wright_fisher_diffusion(double r, std::uint64_t * seed){
      double y;
      if (r > 1e-9) { // TODO citation or reason for this cutoff
        int m = number_of_blocks(r, seed);
        y = sample_beta(1 + m, seed);
      } else {
        y = r / 2;
        std::normal_distribution<double> distribution(0.0, sqrt(r * y * (1 - y)));
        auto sample = prn(seed);
        y = std::abs(sample);
      }
      return y;
    }

    //! Evaluate Full Spherical Brownian Motion
    //!
    //! Evaluates spherical Brownian motion (combined effect of many small scatterings) for an incident proton in the current material
    //! \param distance Distance to evaluate scatter over, i.e. distance incident particle has travelled in current step
    //! \param energy Energy of incident particle
    //! \param direction_in Initial direction of incident particle
    //! \param moliere_transformed_precomp Pre-calculated coeffcient for current material
    //! \param seed Current seed for RNG
    inline std::pair<double, double> spherical_bm(double distance, double energy, std::vector<double> direction_in, double moliere_sd_sq, uint64_t * seed){
      std::vector<double> z, u, w;
      u.resize(3);
      w.resize(3);

      z = direction_in;

      auto y = wright_fisher_diffusion(moliere_sd_sq, seed);
      auto theta = 2.0 * PI * prn(seed);

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
        // might not need the denominator if w unity vector
        return{cos(direction_out_1), direction_out_2};
    }

    //! Evaluate a single elastic scatter event
    //!
    //! Evaluate an elastic scatter event for a proton against the given nuclide at the given energy
    //! \param nuclide The Nuclide to collide against
    //! \param e Energy of incident particle
    //! \param seed Current seed for RNG
    inline double rutherford_elastic_scatter(const Nuclide & nuclide, double e, std::uint64_t * seed){
      auto alpha = nuclide.proton_el_xsec.sample(e*eVToMeV, prn(seed));
      return cos(alpha); //TODO URGENT cm to lab??
    }

    //! Compute separation energies for incident particle
    //!
    //! Computes S_a as in Section 6.2.3.2 on p. 137 [4]
    //! \cite [4] https://doi.org/10.2172/1425114
    //! \param a Atomic Weight of Nuclide
    //! \param z Atomic Number of Nuclide
    //! \return S_a
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

    //! Sample nonelastic collision
    //!
    //! Samples an inelastic collision against a specific nuclide, N
    //! References:
    //! [1] https://doi.org/10.1088/1361-6560/ae5586
    //! [4] https://doi.org/10.2172/1425114
    //! \param nuclide The nuclide to collide with, N
    //! \param e Energy of incident proton, will be updated
    //! \param alpha Polar angle of incident proton, will be updated
    //! \param u A uniform random variate in [0,1]
    //! \param u2 A uniform random variate in [0,1]
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

    //! Evaluate an inelastic collision event
    //!
    //! Evaluates a single inelastic collision with a single Nuclide for a proton at energy e.
    //! 
    //! \param nuclide The Nuclide
    //! \param e Energy of the incident proton in eV
    //! \return A pair, containing the updated energy in eV and the polar scattering angle
    inline std::pair<double, double> non_elastic_scatter(const Nuclide & nuclide, double e, std::uint64_t * seed){
      double alpha;
      e = e*eVToMeV; // e passed by value so working with a COPY below
      sample_nonelastic_collision(nuclide, e, alpha, prn(seed), prn(seed));
      return {e*MeVToeV, alpha};
    }
  };
};
#endif