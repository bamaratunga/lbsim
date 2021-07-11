#ifndef __DEFINITIONS_HPP__
#define __DEFINITIONS_HPP__

#include <cmath>
#include <string>

/******************
 *   6   2   5
 *     \ | /
 *   3 - 0 - 1
 *     / | \
 *   7   4   8
******************/
// D2Q9 LATTICE CONSTANTS
static const int N_DIM = 2;
static const int N_DIRECTIONS = 9;
// Index of lattice velocity in each opposite direction
static const int OPP[N_DIRECTIONS] = {0,   3,   4,   1,   2,   7,   8,   5,   6};
static const int LATTICE_VELOCITIES[N_DIRECTIONS][N_DIM] =  {{ 0, 0},
                                                             { 1, 0},
                                                             { 0, 1},
                                                             {-1, 0},
                                                             { 0,-1},
                                                             { 1, 1},
                                                             {-1, 1},
                                                             {-1,-1},
                                                             { 1,-1}};
static const double LATTICE_WEIGHTS[N_DIRECTIONS] = {4.0/9.0,
                                                     1.0/9.0,
                                                     1.0/9.0,
                                                     1.0/9.0,
                                                     1.0/9.0,
                                                     1.0/36.0,
                                                     1.0/36.0,
                                                     1.0/36.0,
                                                     1.0/36.0};
static const double CS = 1.0/sqrt(3.0);

typedef struct Node {
    double dir[N_DIRECTIONS];
} Node;

typedef struct Vec {
    double comp[N_DIM];
} Vec;

typedef struct Inputs {
    size_t x_dim;
    size_t y_dim;
    size_t timesteps;
    size_t plot_interval;
    bool output_results;

    double reynolds;
    double wall_vel;

    std::string case_name;
    std::string dict_name;
} Inputs;

#endif // __DEFINITIONS_HPP__
