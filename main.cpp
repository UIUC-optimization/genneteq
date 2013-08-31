/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-06-16                                                          */
/* File: main.cpp                                                            */
/* Description:                                                              */
/*   Contains implementation details for the network simplex driver.         */
/* Usage: ./genneteq <args>                                                  */
/*****************************************************************************/
#include "main.h"
#include "matrix.h"
#include "graph.h"
#include "variables.h"
#include "trees.h"
#include "genneteq.h"

#include "getopt.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>

int phaseIpivot;
int phaseIIpivot;

int main(int argc, char **argv) 
{
    Config *conf = new Config();
    parseArgs(conf, argc, argv);

    Graph *g = new Graph();
    (*g).initialize(conf);
    //(*g).print();

    double start = (double)clock();
    GenNetEq *solver = new GenNetEq();

    // Phase 0
    (*solver).initialize(g, conf->initType);
    delete g;    // No longer needed now
    delete conf; // No longer needed now
    double pre = (double)clock();
    printf("GenNetEq finished Phase 0 in %0.3fs\n-----\n",
        (pre - start) / CLOCKS_PER_SEC);

    // Phase I
    int pI_iters = (*solver).solve(USE_PRESOLVE);
    double mid = (double)clock();
    printf("GenNetEq finished Phase I in %0.3fs\n-----\n", 
        (mid - pre) / CLOCKS_PER_SEC);

    // Phase II
    int pII_iters = (*solver).solve(USE_SOLVE);
    double end = (double)clock();
    printf("GenNetEq finished Phase II in %0.3fs\n", 
        (end - mid) / CLOCKS_PER_SEC);

    int totalIters = pI_iters + pII_iters;
    double totalTime = (end - start) / CLOCKS_PER_SEC;
    double pivotTime = (end - pre) / CLOCKS_PER_SEC;
    
    printf("GNE total time: %0.3fs\n", totalTime);
    printf("Pivoting stats: %d, %0.3fs, %0.10f\n", totalIters, pivotTime, 
                                            pivotTime / ((double) totalIters));

    // Clean-up
    delete solver;

    return 0;
}

void parseArgs(Config *conf, int argc, char **argv)
{
    phaseIpivot = LVW; // LVA;
    phaseIIpivot = LVW; // LVA;
    conf->initType = INFEASIBLE;
    int opt;

    // Add any flags to the third argument of getopt; a colon indicates the 
    // option takes an argument, which can be accessed in the variable optarg.
    while ((opt = getopt(argc, argv, "t:p:1:2:?")) != -1) {
        switch(opt) {
        case 't':
            sscanf(optarg, "%d", &(conf->initType));
            break;
        case 'p':
            sscanf(optarg, "%d", &(conf->debug));
            break;
        case '1':
            sscanf(optarg, "%d", &phaseIpivot);
            break;
        case '2':
            sscanf(optarg, "%d", &phaseIIpivot);
            break;
        case '?':
        default:
            printUsageAndExit();
        }
    }

    // optind is the value of the first non-option element of argc.  If it is
    // less than the total number of arguments, then we take in a filename.
    // If there is more than one, print the usage statement and die.
    if (optind == argc - 1) {
        conf->inFile = argv[optind];
    } else { 
        printUsageAndExit();
    }
}

void printUsageAndExit()
{
    printf("Usage: ./genneteq <filename> [-p <debug level>] [-t <int>]\n");
    printf("\t-t  Initial solution type:\n");
    printf("\t    Self loops - %d (default)\n", SELFLOOPS);
    printf("\t    Equal flow self loops - %d\n", EQSELFLOOPS);
    printf("\t    Infeasible initial solution - %d\n", INFEASIBLE);
    printf("\t-1,2  Pivot rule for phase I, II:\n");
    printf("\t    Largest violation - %d\n", LV);
    printf("\t    Largest weighted violation (default) - %d\n", LVW);
    printf("\t    Largest violation (arcs first) - %d\n", LVA);
    printf("\t    Largest violation (eqfs first) - %d\n", LVE);
    printf("\t    First violation - %d\n", FV);
    printf("\t    Steepest Descent - %d\n", SD);
    printf("\t    Random violation - %d\n", RV);
    printf("\t    Random weighted violation - %d\n", RWV);
    printf("\t    Random large violation - %d\n", RLV);
    printf("\t    Random large violation (weighted) - %d\n", RLVW);
    printf("\t    Candidate List - %d\n", CL);
    printf("\t    Candidate List (weighted) - %d\n", CLW);
    printf("\t    Simulated Annealing Rule - (arcs first) %d\n", SAA);
    printf("\t    Simulated Annealing Rule - (eqfs first) %d\n", SAE);
    printf("\t    Candidate List 2 - %d\n", CL2);
    printf("\t    Candidate List (weighted) 2 - %d\n", CLW2);
    printf("\t    Quick Select - %d\n", QS);
    printf("\t    Quick Select (weighted) - %d\n", QSW);
    exit(1);
}

