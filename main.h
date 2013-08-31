/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-06-16                                                          */
/* File: main.h                                                              */
/* Description:                                                              */
/*   Contains design details for the network simplex driver.                 */
/*****************************************************************************/
#ifndef MAIN_H
#define MAIN_H

#include <limits>

#define SELFLOOPS 0
#define EQSELFLOOPS 1
#define INFEASIBLE 2

#define LV 0
#define LVW 1
#define LVA 2
#define LVE 3
#define FV 4
#define SD 5
#define RV 6
#define RWV 7
#define RLV 8
#define RLVW 9
#define CL 10
#define CLW 11
#define SAA 12
#define SAE 13
#define CL2 14
#define CLW2 15
#define QS 16
#define QSW 17

const double ROUNDOFF_TOLERANCE = 
//  1.0e-8;     // 0.00000001           Old value that worked
//  1.0e-10;    // 0.0000000001         New value to try
//  1.0e-12;    // 0.000000000001       New value that worked for most cases
    1.0e-13;    // 0.0000000000001      Further refinement, worked for bad.raw
//  1.0e-14;    // 0.00000000000001    
//  1.0e-15;    // 0.000000000000001    Failed for bad.raw - too tight
//  std::numeric_limits<double>::epsilon(); // Ideal, but not enough slack

const int MAX_ITERATIONS = std::numeric_limits<int>::max();

const bool USE_PRESOLVE = true;
const bool USE_SOLVE = false;

// Contains program options
struct Config 
{
	int initType;
    int debug;
    const char *inFile;
};

// Main functions
void parseArgs(Config *conf, int argc, char **argv);
void printUsageAndExit();

#endif // MAIN_H

