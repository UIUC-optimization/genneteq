/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-06-16                                                          */
/* File: genneteq.cpp                                                        */
/* Description:                                                              */
/*   Contains implementation details for the GenNetEq algorithm.             */
/*****************************************************************************/
#include "main.h"
#include "variables.h"
#include "trees.h"
#include "genneteq.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include <vector>
using std::vector;

GenNetEq::GenNetEq()
{
    // Use for simulated annealing only, just set to some random big number.
    lastViolation = 100000;
    baselineViolation = 100000;

	totalPivots = 0;
	arcArcPivots = 0;
	arcEFSPivots = 0;
	EFSArcPivots = 0;
	EFSEFSPivots = 0;
	sameTIPivots = 0;
	sameTIIPivots = 0;
	tItIPivots = 0;
	tIItIIPivots = 0;
	tItIIPivots = 0;
	EFSPivots = 0;
	avgTIsize = 0.0;
	avgTIIsize = 0.0;
}

GenNetEq::~GenNetEq()
{
    clearCandidateList();
    while (!cycleList.empty()) {
        delete cycleList.front();
        cycleList.pop_front();
    }
    for (int i = 0; i < setBtI.size(); ++i) {
        delete setBtI[i];
    }
    for (int i = 0; i < setBtII.size(); ++i) {
        delete setBtII[i];
    }
    for (int varInd = 0; varInd < vars.size(); ++varInd) {
        delete vars[varInd];
    }
}

int GenNetEq::solve(bool usePresolve)
{
    // Initialize statistics and variables
    bool pIwasInfeasible = isInfeasible;
    initializeStats(usePresolve);               // O(1)
    setVariableCosts();

    double start = (double)clock();

    // If we aren't presolving, then the basis we finished Phase I with should 
    // be feasible. If it isn't, we terminate since the problem is infeasible. 
    if (!presolve && pIwasInfeasible) {
        isInfeasible = true;
        printTerminationStatus();
        return loopIter;
    }

    // Compute initial flow (also computes initial cost)
    computeFlows();                             // O(m' + np + g(s))
//    double gain = 0.0;

    // Begin main loop
    while ((!isOptimal) && (!isUnbounded) && 
           ((phaseIiters + phaseIIiters) < MAX_ITERATIONS)) {
//        printf("Iter %6d; Cost: %lf; Gain: %lf\n", loopIter, flowCost, gain);

        if ((loopIter % 1000) == 0) {
            // Recompute flows from scratch to reduce roundoff error 
            // accumulation
            computeFlows(2); // Check and update
            printf("Iter %7d (%7.2lfs): Obj: %lf\n", loopIter, 
                (clock() - start) / CLOCKS_PER_SEC, flowCost);
//        } else {
//            computeFlows(1); // Just check
        }
//        double oldCost = flowCost;

        computePotentials();                        // O(ns + g(s))
        identifyEnteringVariable();                 // O(t_select) = O(m' + np)
        if (eVar != NULL) {
            pivot();                                // O(ns + g(s))
            ++loopIter;
            presolve ? ++phaseIiters : ++phaseIIiters;

        } // else we should be optimal
//        gain = oldCost - flowCost;
    }

    computeFlows();
    printf("Iter %7d (%7.2lfs): Obj: %lf\n", loopIter, 
        (clock() - start) / CLOCKS_PER_SEC, flowCost);

    if (presolve) checkFeasibility();           // O(m' + n + p)
    verifyOptimality();                         // O(ns + g(s) + m' + np)
    printTerminationStatus();

    return loopIter;
}

/*****************************************************************************/
/* Miscellaneous Protected Functions                                         */
/*****************************************************************************/
// Running time: O(1)
void GenNetEq::initializeStats(bool usePresolve)
{
    presolve = usePresolve;

    // Standard stat and progress counters
    isOptimal = false;
    isInfeasible = false;
    isUnbounded = false;

    loopIter = 0;
    infeasLoopIter = 0;

    consecutiveCycles = 0;
    consecutiveDegenPivots = 0;
    totalCyclesDetected = 0;
    totalDegenPivots = 0;

    cyclingDetected = false;
    lastEVarInd = -1;
    itersSinceUpdate = -1;
    // Parameter Settings
    trackLowestVarIndex = false;
    selectRandomBlockingVar = true;
    selectRandomOffendingVar = true;
    useRandomJumpWhenCycling = true;
    numLargestViolVarsToTrack = 5;
    majorIterationLength = 3;

    clMaxSize = numNodes;
    nbArcPivotInd = 0;
    nbEqfPivotInd = 0;
    numArcsToAlwaysCheck = 5;
    numEqfsToAlwaysCheck = 1;

    //if ((presolve && phaseIiters == 0) || phaseIIiters == 0)
    baselineViolation = lastViolation;

    // Violation threshold is also a parameter, used for scaled pivoting
    if (presolve) {
        violationThreshold = (double) numNodes; // Just a guess?
    } else {
        violationThreshold = bigM * ((double) numNodes);
    }
    return;
}

// Running time: O(m' + p)
void GenNetEq::printTerminationStatus()
{
    if (presolve) {
        printf("Finished GenNetEq Phase I in %d iterations\n", loopIter);
    } else {
        printf("Finished GenNetEq Phase II in %d iterations\n", loopIter);
    }

    // Determine termination status
    if (isUnbounded) {
        printf("Problem is unbounded.\n");
    } else if (isInfeasible) {
        printf("Problem is infeasible.\n");
    } else { // isOptimal
        printf("Problem is feasible.\n");
        //computeFlowCost();                  // O(m' + p)
        if (presolve) {
            printf("The initial cost is %lf\n", flowCost);
        } else {
            printf("The optimal cost is %lf\n", flowCost);
        }
        //printf("The optimal basis is:\n");
        //prettyPrintBasis();
        //prettyPrintFlow();
    }
    printf("There were %d degenerate pivots\n", totalDegenPivots);
    printf("There were %d cycles detected\n", totalCyclesDetected);

/*
	printf("%d/%d arc-arc pivots (%0.1f\%)\n", arcArcPivots, totalPivots, 
			(arcArcPivots * 100.0) / totalPivots);
	printf("%d/%d arc-efs pivots (%0.1f\%)\n", arcEFSPivots, totalPivots, 
			(arcEFSPivots * 100.0) / totalPivots);
	printf("%d/%d efs-arc pivots (%0.1f\%)\n", EFSArcPivots, totalPivots, 
			(EFSArcPivots * 100.0) / totalPivots);
	printf("%d/%d efs-efs pivots (%0.1f\%)\n\n", EFSEFSPivots, totalPivots, 
			(EFSEFSPivots * 100.0) / totalPivots);
	printf("%d/%d entering vars in a single TI tree (%0.1f\%)\n", sameTIPivots, totalPivots, 
			(sameTIPivots * 100.0) / totalPivots);
	printf("Average size of type I trees: %0.1f\n", (avgTIsize + 0.0)/ sameTIPivots);
	printf("%d/%d entering vars in a single TII tree (%0.1f\%)\n", sameTIIPivots, totalPivots, 
			(sameTIIPivots * 100.0) / totalPivots);
	printf("Average size of type II trees: %0.1f\n", (avgTIIsize + 0.0)/ sameTIIPivots);
	printf("%d/%d entering vars cross two TI trees (%0.1f\%)\n", tItIPivots, totalPivots, 
			(tItIPivots * 100.0) / totalPivots);
	printf("%d/%d entering vars cross two TII trees (%0.1f\%)\n", tIItIIPivots, totalPivots, 
			(tIItIIPivots * 100.0) / totalPivots);
	printf("%d/%d entering vars go between TI and TII trees (%0.1f\%)\n", tItIIPivots, totalPivots, 
			(tItIIPivots * 100.0) / totalPivots);
	printf("%d/%d entering vars are equal flow sets (%0.1f\%)\n", EFSPivots, totalPivots,
			(EFSPivots * 100.0) / totalPivots);
*/

    return;
}

// Running time: O(m' + n + p)
void GenNetEq::checkFeasibility() 
{
    // Determine if the current basis system represents a feasible flow (used 
    // to stop Phase I). By default we assume feasibility and then search for 
    // violations of it. 
    isInfeasible = false;

    // If any of the arcs have ``fake'' costs associated with them (because 
    // the flow is infeasible), the basis isn't feasible.
    for (int i = 0; i < useArtificialCosts.size(); ++i) {
        if (useArtificialCosts[i] != 0.0) {
            isInfeasible = true;
            printf("Arc has associated artificial cost\n");
            return;
        }
    }
    
    // If any basic arc or equal flow set is artificial with non-zero flow or 
    // violates its capacity constraints, problem is infeasible
    for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
        ArcVar *xb = setBarc[bArcInd];
        if (((xb->isArtificial) && (xb->flow > ROUNDOFF_TOLERANCE)) ||
            (xb->flow < xb->LB - ROUNDOFF_TOLERANCE) ||
            (xb->flow > xb->UB + ROUNDOFF_TOLERANCE)) {
            isInfeasible = true;
            printf("Artificial arc in basis with non-zero flow %.16lf\n", 
                    xb->flow);
            return;
        }
    }
    for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
        EqFlowVar *xb = setBeqf[bEqfInd];
        if (((xb->isArtificial) && (xb->flow > ROUNDOFF_TOLERANCE)) ||
            (xb->flow < xb->LB - ROUNDOFF_TOLERANCE) ||
            (xb->flow > xb->UB + ROUNDOFF_TOLERANCE)) {
            isInfeasible = true;
            printf("Artificial equal flow set in basis with non-zero flow "
                   "%.16lf\n", xb->flow);
            return;
        }
    }

    // If any artificial variable is at its upper bounds, problem is infeasible
    for (int nbArcInd = 0; nbArcInd < setNBarc.size(); ++nbArcInd) {
        ArcVar *xnb = setNBarc[nbArcInd];
        if (xnb->isArtificial && xnb->setLoc == U_var) {
            isInfeasible = true;
            printf("Artificial arc at upper bound\n");
            return;
        }
    }
    for (int nbEqfInd = 0; nbEqfInd < setNBeqf.size(); ++nbEqfInd) {
        EqFlowVar *eqfv = setNBeqf[nbEqfInd];
        if (eqfv->isArtificial && eqfv->setLoc == U_var) {
            isInfeasible = true;
            printf("Artificial equal flow set at upper bound\n");
            return;
        }
    }

    return;
}

// Running time: O(m' + p)
void GenNetEq::computeFlowCost()
{
    double newFlowCost = 0.0;
    for (int i = 0; i < vars.size(); ++i) {
        newFlowCost += vars[i]->flow * vars[i]->cost;
    }
    printf("Current Flow Cost is %lf, actual flow cost should be %lf\n", 
            flowCost, newFlowCost);
    flowCost = newFlowCost;
    return;
}

// Running time: O(m' + p)
void GenNetEq::setVariableCosts()
{
    if (presolve) {
        // Total flow on artificial variables
        if (initialType == INFEASIBLE) {
            for (int i = 0; i < vars.size(); ++i) {
                vars[i]->cost = vars[i]->isArtificial ? 1.0 : useArtificialCosts[i]; 
            }
        } else { 
            for (int i = 0; i < vars.size(); ++i) {
                vars[i]->cost = vars[i]->isArtificial; 
            }
        }
    } else {
        for (int i = 0; i < vars.size(); ++i) {
            vars[i]->cost = vars[i]->actualCost;
        }
    }
    return;
}

