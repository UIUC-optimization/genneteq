/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-06-16                                                          */
/* File: gnePotential.cpp                                                    */
/* Description:                                                              */
/*   Contains implementation details for the potential functions used by the */
/*   GenNetEq algorithm.                                                     */
/*****************************************************************************/
#include "main.h"
#include "variables.h"
#include "trees.h"
#include "genneteq.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <cstdio>
#include <cstdlib>

#include <vector>
using std::vector;

/*****************************************************************************/
/* Potentials Protected Functions                                            */
/*****************************************************************************/
// Running time: O(ns + g(s))
void GenNetEq::computePotentials()
{
    // Initialize potentials in terms of root potentials [O(n)]
    for (int tIind = 0; tIind < setBtI.size(); ++tIind) {
        computePotentialsInTermsOfRoot(setBtI[tIind]);
        finalizePotentialsTypeI(setBtI[tIind]);
    }
    for (int tIIind = 0; tIIind < setBtII.size(); ++tIIind) {
        computePotentialsInTermsOfRoot(setBtII[tIIind]);
    }

    if (setBeqf.size() > 0) {   
        // Now that we have the potentials with variables, we can set up the 
        // system of equations to solve for variables [O(ns + g(s))]
        computePotentialsWithEqFlowSets();  
    }
    // Otherwise all nodes belong to type I trees, and hence all have 
    // their potentials already computed and moved to the appropriate field. 

    // DEPRECATED:
//    } else {
//        // All potentials are already computed, but we just need to move them 
//        // over to the appropriate potentials vector [O(n)]
//        for (int i = 0; i < numNodes; ++i) {
//            nodes[i].potential = nodes[i].potentialWithEq.a0;
//            // I don't remember writing this...
//            //if (fabs(nodes[i].potentialWithEq.a0) > ROUNDOFF_TOLERANCE)
//            //    nodes[i].potential = ((*potentialsWithEqs)[i])->a0;
//            //else
//            //    nodes[i].potential = 0.0;
//        }
//    }
    return;
}

// Running time: O(|Ti|)
void GenNetEq::computePotentialsInTermsOfRoot(Tree *tX)
{
    int treeThetaID; 
    if ((*tX).getTreeType() == TYPE_I) {
        treeThetaID = THETA;
    } else { // (*tX).getTreeType() == TYPE_II
        treeThetaID = (*tX).getTreeID();
    }

    (*tX).resetRootToLeafTraversal();
    int h = (*tX).getNextThread(); // This is the root

    ThetaEquation &rootNpEq = nodes[h].potentialWithEq;
    // pi_{h} = theta
    rootNpEq.a0 = 0.0;
    rootNpEq.a1 = 1.0;     // Root's potential equals 1.0 * theta[treeThetaID]
    rootNpEq.thetaID = treeThetaID;

    // Now run through tree in forward thread order, setting potentials as 
    // nodes are visited [O(|Ti|)]
    int j = (*tX).getNextThread();
    while (j != h) {
        int i = nodes[j].pred;
        ArcVar *predArc = nodes[j].predArc;

        ThetaEquation &npEqI = nodes[i].potentialWithEq;
        ThetaEquation &npEqJ = nodes[j].potentialWithEq;

        double c = predArc->cost;
        double mu = predArc->mult;

        if (predArc->tail == i) {
            // pi(j) = (pi(i) - c_ij) / mu_ij
            npEqJ.a0 = (npEqI.a0 - c) / mu;
            npEqJ.a1 = (npEqI.a1) / mu;
            npEqJ.thetaID = treeThetaID;
        } else { // (predArc->tail == j) 
            // pi(j) = mu_ji pi(i) + c_ji
            npEqJ.a0 = (mu * (npEqI.a0)) + c;
            npEqJ.a1 = mu * (npEqI.a1);
            npEqJ.thetaID = treeThetaID;
        }

        j = (*tX).getNextThread();
    }
    return;
}

// Running time: O(|Ti|)
void GenNetEq::finalizePotentialsTypeI(TreeTypeI *tI)
{
    // Compute value of theta
    ArcVar *extraArc = (*tI).getExtraArc();
    int alpha = extraArc->tail;
    int beta = extraArc->head;
    double cAB = extraArc->cost;
    double muAB = extraArc->mult;
    ThetaEquation &npEqAlpha = nodes[alpha].potentialWithEq; 
    ThetaEquation &npEqBeta = nodes[beta].potentialWithEq;
    double thetaValue = (cAB - (npEqAlpha.a0) + muAB * (npEqBeta.a0)) 
                            / ((npEqAlpha.a1) - muAB * (npEqBeta.a1));

    // Now compute actual potentials for all nodes in the tree [O(|Ti|)]
    vector<int> *tNodes = (*tI).getNodesInTree();
    for (int tnInd = 0; tnInd < (*tNodes).size(); ++tnInd) {
        int i = (*tNodes)[tnInd];
        ThetaEquation &npEqI = nodes[i].potentialWithEq;
        nodes[i].potential = npEqI.a0 + npEqI.a1 * thetaValue;
        npEqI.a0 = 0.0;
        npEqI.a1 = 0.0;
        npEqI.thetaID = NULL_THETA;
    }

    return;
}

// Running time: O(ns + g(s)) 
// g(s) := Time for GNU GSL to solve a system of s equations in s unknowns
void GenNetEq::computePotentialsWithEqFlowSets()
{
    int temp;
    int s = setBeqf.size();
    gsl_vector *x = gsl_vector_alloc(s);
    gsl_vector *rhs = gsl_vector_alloc(s);
    gsl_matrix *mat = gsl_matrix_alloc(s, s);
    gsl_permutation *p = gsl_permutation_alloc(s);

    // DEBUG:
//    printBasis();
//    printTrees();
//    printPotentialsWithEqs();
//    fflush(stdout);

    // Set up values [O(ns + s^2) <= O(ns)]
    for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
        EqFlowVar *eqfv = setBeqf[bEqfInd];
        // Calculate values for this equation line
        double rhsVal = eqfv->cost;
        vector<double> *rowVals = new vector<double>(s, 0.0);
        for (int i = 0; i < numNodes; ++i) {
            double dri = eqfv->nodeVals[i];
            if (nodes[i].memberOfTreeType == TYPE_I) {
                rhsVal -= dri * nodes[i].potential;
            } else { // nodes[i].memberOfTreeType == TYPE_II
                ThetaEquation &npEq = nodes[i].potentialWithEq;
                rhsVal -= dri * (npEq.a0);
                (*rowVals)[npEq.thetaID] += dri * (npEq.a1); 
            }
        }
        // Update values in gsl objects
        gsl_vector_set(rhs, bEqfInd, rhsVal);
        for (int tIIind = 0; tIIind < s; ++tIIind) {
            gsl_matrix_set(mat, bEqfInd, tIIind, (*rowVals)[tIIind]);
        }
        delete rowVals;
    }

    // DEBUG:
//    printf("Getting ready to solve system of equations\n");
//    printGSLmatrix(mat, s, s);
//    fflush(stdout);

    // Solve the system [O(g(s))]
    gsl_linalg_LU_decomp(mat, p, &temp);
    gsl_linalg_LU_solve(mat, p, rhs, x);

    // Now values in x correspond to the node potentials for the root nodes 
    // of the type II trees, so we can use them to calculate numeric 
    // potentials for all these nodes [O(n)]
    // NOTE: This could be made tighter by iterating over the type II trees 
    // and the nodes in those trees, but I think memory accesses this way 
    // might be more favorable for speed
    for (int i = 0; i < numNodes; ++i) {
        if (nodes[i].memberOfTreeType == TYPE_I) { // Already solved for them
            continue;
        }
        ThetaEquation &npEq = nodes[i].potentialWithEq;
        nodes[i].potential = npEq.a0 + npEq.a1 * gsl_vector_get(x,npEq.thetaID);
        npEq.a0 = 0.0;
        npEq.a1 = 0.0;
        npEq.thetaID = NULL_THETA;
    }

    // Clean up gsl components
    gsl_permutation_free(p);
    gsl_matrix_free(mat);
    gsl_vector_free(rhs);
    gsl_vector_free(x);

    return;
}

// Running time: O(m' + np)
void GenNetEq::computeReducedCosts(bool includeBasicVars) 
{
    // By default assume that the basis is optimal until we find a violation
    isOptimal = true;
    offendingVars.clear();

    // Compute reduced costs and check optimality conditions for arcs at 
    // lower and upper bounds [O(m')] 
    for (int nbArcInd = 0; nbArcInd < setNBarc.size(); ++nbArcInd) {
        ArcVar *nba = setNBarc[nbArcInd];
        double rc = computeArcRedCost(nba);
        if (((nba->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((nba->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            offendingVars.push_back(nba);
        }
    }

    // Compute reduced costs and check optimality conditions for equal flow 
    // sets at lower and upper bounds [O(np)] 
    for (int nbEqfInd = 0; nbEqfInd < setNBeqf.size(); ++nbEqfInd) {
        EqFlowVar *nbeq = setNBeqf[nbEqfInd];
        double rc = computeEqfRedCost(nbeq);
        if (((nbeq->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((nbeq->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            offendingVars.push_back(nbeq);
        }
    }

    if (includeBasicVars) {
        // Compute reduced costs for basic arcs and equal flow sets [O(ns)]
        for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
            computeArcRedCost(setBarc[bArcInd]);
        }
        for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
            computeEqfRedCost(setBeqf[bEqfInd]);
        }
    }

    return;
}

// Running time: O(m' + np + g(s))
void GenNetEq::verifyOptimality()
{
    // First make sure that potentials and reduced costs are up to date
    computePotentials();
    computeReducedCosts(true);

    if (!offendingVars.empty()) {
        printf("%d non-basic variables violate their optimality conditions\n", 
            offendingVars.size());
    }

    // Check optimality conditions for basic vars
    for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
        if (fabs(setBarc[bArcInd]->redCost) > ROUNDOFF_TOLERANCE) {
            printf("Basic var has non-zero red. cost: %.16lf\n", 
                setBarc[bArcInd]->redCost);
        }
    }
    for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
        if (fabs(setBeqf[bEqfInd]->redCost) > ROUNDOFF_TOLERANCE) {
            printf("Basic var has non-zero red. cost: %.16lf\n", 
                setBeqf[bEqfInd]->redCost);
        }
    }

    return;
}

