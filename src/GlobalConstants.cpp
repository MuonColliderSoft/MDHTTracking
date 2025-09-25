#include "include/GlobalConstants.h"

// Fundamental constants
const double Pi = 3.14159265;
// Speed of light: 299.792458 mm/ns
const double SpeedOfLight = 299.792458;

// Debugging flags
bool Special = false;
bool verbose = false;

// PARTICLE MASSES in GeV/c2
const unsigned nMasses = 5;
double Mass[] = {0., 0.000511, 0.139570, 0.493677, 0.938272};

// Training phase (global)
int TrainingPhase = 0;
