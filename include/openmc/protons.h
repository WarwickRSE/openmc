#ifndef __proton_physics__
#define __proton_physics__

#include <random>
//TODO - what RNG to use?? NOT a file local static one pls!!

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



#endif