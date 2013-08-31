/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-06-16                                                          */
/* File: gneFlow.cpp                                                         */
/* Description:                                                              */
/*   Contains implementation details for the flow functions used by the      */
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
/* Flows Protected Functions                                                 */
/*****************************************************************************/
// Running Time:
//   Using Equal Flow Sets:         O(m' + np + g(s))
//   Not Using Equal Flow Sets:     O(m' + np)
void GenNetEq::computeFlows(int updateOption)
{
    flowCost = 0.0;
    if (setBeqf.size() > 0) {
        // Compute flows for all arcs in type II trees and for the equal flow 
        // sets, plus update supply for all nodes [O(m' + np + g(s))]
        setSupplyWithThetas();
        computeBoundedFlows(true);
        computeFlowsWithEqFlowSets();
    } else {
        // No type II trees or equal flow sets [O(m' + np)]
        setSupplyWithEqs();
        computeBoundedFlows(false);
    }

    // At this point, we have numerical flow values for all basic equal flow 
    // sets and type II tree arcs, and supplyWithEqs has been set 
    // appropriately with numerical values for all nodes in type I trees and 
    // reset to 0 for all nodes in type II trees. Now, apply the procedure to 
    // compute flow on the type I trees. [O(n)]
    for (int tIind = 0; tIind < setBtI.size(); ++tIind) {
        computeFlowsTypeI(setBtI[tIind]);
    }

    setFlowValues(updateOption);
    return;
}

// Running time: O(m' + p)
void GenNetEq::setFlowValues(int updateOption)
{
    // Now we assume that the new flow values for the basic variables are 
    // contained in the deltaFlow fields, so we can copy the flow values over 
    // to the standard flow field or apply a check as desired. [O(m' + p)]
    if (updateOption == 0) { // Apply standard recompute
        for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
            ArcVar *av = setBarc[bArcInd];
            av->flow = av->deltaFlow;
            av->deltaFlow = 0.0;
            flowCost += av->flow * av->cost;
        }
        for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
            EqFlowVar *eqfv = setBeqf[bEqfInd];
            eqfv->flow = eqfv->deltaFlow;
            eqfv->deltaFlow = 0.0;
            flowCost += eqfv->flow * eqfv->cost;
        }
    } else { // Check for accumulated roundoff error
        int numErrors = 0;
        int curFlowInfeasibilities = 0;
        int recompFlowInfeasibilities = 0;
        double diffMax = 0.0;

        // First check the basic arcs
        for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
            ArcVar *av = setBarc[bArcInd];
            if (fabs(av->deltaFlow - av->flow) > 1.0e-6) {
                ++numErrors;
                if (diffMax < fabs(av->deltaFlow - av->flow)) {
                    diffMax = fabs(av->deltaFlow - av->flow);
                }
                //printf("Recomp Flow: %.16lf != %.16lf :Cur Flow (%.16lf)\n", 
                //     av->deltaFlow, av->flow, av->deltaFlow - av->flow);
            }
            if (!checkVarFlow(av, false)) {
                // Check current flow feasibility
                ++curFlowInfeasibilities;
            }
            if (!checkVarFlow(av, true)) {
                // Check recomputed flow feasibility
                ++recompFlowInfeasibilities;
            }
            if (updateOption > 1) { // Actually update in addition to checking
                av->flow = av->deltaFlow; 
            }
            flowCost += av->flow * av->cost;
            av->deltaFlow = 0.0;
        }

        // Now check the basic equal flow sets
        for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
            EqFlowVar *eqfv = setBeqf[bEqfInd];
            if (fabs(eqfv->deltaFlow - eqfv->flow) > 1.0e-6) {
                ++numErrors;
                if (diffMax < fabs(eqfv->deltaFlow - eqfv->flow)) {
                    diffMax = fabs(eqfv->deltaFlow - eqfv->flow);
                }
                //printf("Recomp Flow: %.16lf != %.16lf :Cur Flow (%.16lf)\n", 
                //     eqfv->deltaFlow, eqfv->flow, eqfv->deltaFlow - eqfv->flow);
            }
            if (!checkVarFlow(eqfv, false)) {
                // Check current flow feasibility
                ++curFlowInfeasibilities;
            }
            if (!checkVarFlow(eqfv, true)) {
                // Check recomputed flow feasibility
                ++recompFlowInfeasibilities;
            }
            if (updateOption > 1) { // Actually update in addition to checking
                eqfv->flow = eqfv->deltaFlow; 
            }
            flowCost += eqfv->flow * eqfv->cost;
            eqfv->deltaFlow = 0.0;
        }

        // Now print out the error report
        if (numErrors+curFlowInfeasibilities+recompFlowInfeasibilities > 0) {
            printf("Iteration %d\n", loopIter);
            printf("%d variables with significant (> 1.0e-6) accumulated " 
                   "errors in flow\n", numErrors); 
            printf("Max error is %.16lf\n", diffMax);
            printf("%d current flow values violate feasibility\n", 
                    curFlowInfeasibilities);
            printf("%d recomputed flow values violate feasibility\n", 
                    recompFlowInfeasibilities);
        }
    }

    return;
}

/*****************************************************************************/
/* Functions for computing flow when pivoting                                */
/*****************************************************************************/
// Running time: O(|Tu| + |Tv|)
void GenNetEq::computeFlowsPivotingTI(ArcVar *eav, TreeTypeI *tU, 
                                                   TreeTypeI *tV)
{
    // Set up supply modification for given arc. NOTE: This is essentially 
    // what gets done had we called setSupplyWithEqs(eav).
    if (eav->setLoc == L_var) { // eav at LB
        if (eav->tail == eav->head) {
            nodes[eav->tail].supplyWithEq.a0 = -1.0 * (1.0 - eav->mult);
        } else {
            nodes[eav->tail].supplyWithEq.a0 = -1.0;
            nodes[eav->head].supplyWithEq.a0 = eav->mult;
        }
    } else {
        if (eav->tail == eav->head) {
            nodes[eav->tail].supplyWithEq.a0 = (1.0 - eav->mult);
        } else {
            nodes[eav->tail].supplyWithEq.a0 = 1.0;
            nodes[eav->head].supplyWithEq.a0 = -1.0 * eav->mult;
        }
    }

    // At this point, supplyWithEqs has been set appropriately with numerical 
    // values for all nodes in the appropriate trees. Now, apply the procedure 
    // to compute flow on the type I trees. [O(|Tu| + |Tv|)]
    computeFlowsTypeI(tU);
    if (tU->getTreeID() != tV->getTreeID()) {
        computeFlowsTypeI(tV);
    }

    return;
}

// Running time: O(ns + g(s))
void GenNetEq::computeFlowsPivotingTII(Variable *var)
{
    // Compute flows for all arcs in type II trees and for the equal flow sets, 
    // plus update supply for all nodes [O(ns + g(s))]
    if (setBeqf.size() > 0) {
        setSupplyWithThetas(var);
        computeFlowsWithEqFlowSets();
    } else {
        setSupplyWithEqs(var);
    }

    // At this point, we have numerical flow values for all basic equal flow 
    // sets and type II tree arcs, and supplyWithEqs has been set 
    // appropriately with numerical values for all nodes in type I trees and 
    // reset to 0 for all nodes in type II trees. Now, apply the procedure to 
    // compute flow on the type I trees. [O(n)]
    for (int tIind = 0; tIind < setBtI.size(); ++tIind) {
        computeFlowsTypeI(setBtI[tIind]);
    }

    return;
}

// Running time: O(ns + g(s))
void GenNetEq::computeFlowsWithEqFlowSets()
{
    // Implicitly set up variable flows on basic equal flow sets, and update 
    // variable supplies at all nodes based on d_r(i) [O(ns)]
    for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
        EqFlowVar *eqfv = setBeqf[bEqfInd];
        for (int i = 0; i < numNodes; ++i) {
            nodes[i].supplyWithThetas[bEqfInd + 1] -= eqfv->nodeVals[i]; 
        }
    }

    // Initialize gsl components for solving system O(s^2) <= O(ns)
    int temp;
    int s = setBeqf.size();
    gsl_vector *x = gsl_vector_alloc(s);
    gsl_vector *rhs = gsl_vector_alloc(s);
    gsl_matrix *mat = gsl_matrix_alloc(s, s);
    gsl_permutation *p = gsl_permutation_alloc(s);

    // For each type II tree, compute flows in terms of theta_r's [O(ns)]
    for (int tIIind = 0; tIIind < setBtII.size(); ++tIIind) {
        int h = computeFlowsTypeII(setBtII[tIIind]);
        vector<double> &rootSupply = nodes[h].supplyWithThetas;
        gsl_vector_set(rhs, tIIind, (-1.0) * (rootSupply[0]));
        for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
            gsl_matrix_set(mat, tIIind, bEqfInd, rootSupply[bEqfInd + 1]);
        }
    }

    // Solve the system [O(g(s))]
    gsl_linalg_LU_decomp(mat, p, &temp);
    gsl_linalg_LU_solve(mat, p, rhs, x);
  
    // Now the variables are solved for, so we can move all information over 
    // to flowsWithEqs and supplyWithEqs [O(ns)]
    computeNumericFlowsAndSupply(x);

    // Clean-up
    gsl_permutation_free(p);
    gsl_matrix_free(mat);
    gsl_vector_free(rhs);
    gsl_vector_free(x);

    return;
}

// Running time: O(ns)
void GenNetEq::computeNumericFlowsAndSupply(gsl_vector *x)
{
    // Now values in x correspond to the flow values for the basic equal flow 
    // sets, so we can use them to calculate flow values for all arcs in the 
    // type II trees [O(ns)]
    for (int tIIind = 0; tIIind < setBtII.size(); ++tIIind) {
        vector<ArcVar *> *tIIarcs = setBtII[tIIind]->getArcsInTree();
        for (int taInd = 0; taInd < (*tIIarcs).size(); ++taInd) {
            ArcVar *av = (*tIIarcs)[taInd];
            vector<double> &avFlow = av->flowWithThetas;
            double actualFlow = avFlow[0];
            for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
                actualFlow += avFlow[bEqfInd+1] * gsl_vector_get(x, bEqfInd);
            }
            av->deltaFlow = actualFlow;
            av->flowWithThetas.clear();
            checkVarForBlocking(av); // NEW
        }
    } 

    // Now compute the flows for the actual equal flow sets [O(s)]
    for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
        EqFlowVar *eqfv = setBeqf[bEqfInd];
        eqfv->deltaFlow = gsl_vector_get(x, bEqfInd);
        checkVarForBlocking(eqfv); // NEW
    }

    // Now calculate the actual supply values at all nodes [O(ns)]
    for (int i = 0; i < numNodes; ++i) {
        vector<double> &varNodeSupply = nodes[i].supplyWithThetas;
        double actualSupply = varNodeSupply[0];
        for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
            actualSupply += 
                varNodeSupply[bEqfInd+1] * gsl_vector_get(x, bEqfInd);
        }
        nodes[i].supplyWithEq.a0 = actualSupply;
        if (nodes[i].memberOfTreeType == TYPE_II) {
            // NOTE: The actual supply should be 0 already, but doing this 
            // explicitly will help reduce accumulation of roundoff errors
            nodes[i].supplyWithEq.a0 = 0.0;
        }
        nodes[i].supplyWithEq.a1 = 0.0;
        nodes[i].supplyWithEq.thetaID = NULL_THETA;
        nodes[i].supplyWithThetas.clear();
    }

    return;
}

/*****************************************************************************/
/* Used for standard flow computation (not pivoting)                         */
/*****************************************************************************/
// Running time: O(m' + np)
void GenNetEq::computeBoundedFlows(bool withEqFlowSets)
{
    if (withEqFlowSets) {
        // Handle the non-basic arcs
        for (int nbArcInd = 0; nbArcInd < setNBarc.size(); ++nbArcInd) {
            ArcVar *nba = setNBarc[nbArcInd];
            if (nba->setLoc == L_var) {
                nba->flow = nba->LB;
            } else { // nba->setLoc == U_var
                nba->flow = nba->UB;
                flowCost += nba->flow * nba->cost;
                nodes[nba->tail].supplyWithThetas[0] -= nba->flow;
                nodes[nba->head].supplyWithThetas[0] += nba->flow * nba->mult;
            }
        }
        // Handle the non-basic equal flow sets
        for (int nbEqfInd = 0; nbEqfInd < setNBeqf.size(); ++nbEqfInd) {
            EqFlowVar *nbeq = setNBeqf[nbEqfInd];
            if (nbeq->setLoc == L_var) {
                nbeq->flow = nbeq->LB;
            } else { // nbeq->setLoc == U_var
                nbeq->flow = nbeq->UB;
                flowCost += nbeq->flow * nbeq->cost;
                for (int i = 0; i < numNodes; ++i) {
                    nodes[i].supplyWithThetas[0] -= 
                                                nbeq->flow * nbeq->nodeVals[i]; 
                }
            }
        }
    } else {
        // Handle the non-basic arcs
        for (int nbArcInd = 0; nbArcInd < setNBarc.size(); ++nbArcInd) {
            ArcVar *nba = setNBarc[nbArcInd];
            if (nba->setLoc == L_var) {
                nba->flow = nba->LB;
            } else { // nba->setLoc == U_var
                nba->flow = nba->UB;
                flowCost += nba->flow * nba->cost;
                nodes[nba->tail].supplyWithEq.a0 -= nba->flow;
                nodes[nba->head].supplyWithEq.a0 += nba->flow * nba->mult;
            }
        }
        // Handle the non-basic equal flow sets
        for (int nbEqfInd = 0; nbEqfInd < setNBeqf.size(); ++nbEqfInd) {
            EqFlowVar *nbeq = setNBeqf[nbEqfInd];
            if (nbeq->setLoc == L_var) {
                nbeq->flow = nbeq->LB;
            } else { // nbeq->setLoc == U_var
                nbeq->flow = nbeq->UB;
                flowCost += nbeq->flow * nbeq->cost;
                for (int i = 0; i < numNodes; ++i) {
                    nodes[i].supplyWithEq.a0 -= nbeq->flow * nbeq->nodeVals[i];
                }
            }
        }
    }
    return;
}

/*****************************************************************************/
/* Set up supply values                                                      */
/*****************************************************************************/
// Running time: O(ns)
void GenNetEq::setSupplyWithThetas()
{
    int s = setBeqf.size();
    // Set up variable supplies at nodes based on supply vector [O(ns)]
    for (int i = 0; i < numNodes; ++i) {
        nodes[i].supplyWithThetas.resize(s + 1, 0.0);
        nodes[i].supplyWithThetas[0] = nodes[i].supply;
    }
    return;
}

// Running time: O(n)
void GenNetEq::setSupplyWithEqs()
{
    for (int i = 0; i < numNodes; ++i) {
        nodes[i].supplyWithEq.a0 = nodes[i].supply;
        nodes[i].supplyWithEq.a1 = 0.0;
        nodes[i].supplyWithEq.thetaID = NULL_THETA;
    }
    return;
}

// Running time: O(ns)
void GenNetEq::setSupplyWithThetas(Variable *var)
{
    int s = setBeqf.size();
    if (var->varType == ARC_VAR) {
        ArcVar *eav = dynamic_cast<ArcVar *>(var);
        // Set up variable supplies of 0 at all nodes [O(ns)]
        for (int i = 0; i < numNodes; ++i) {
            nodes[i].supplyWithThetas.resize(s + 1, 0.0);
        }
        // Now update the supplies at the given arc
        if (eav->setLoc == L_var) { // eav at LB
            if (eav->tail == eav->head) {
                nodes[eav->tail].supplyWithThetas[0] = -1.0 * (1.0 - eav->mult);
            } else {
                nodes[eav->tail].supplyWithThetas[0] = -1.0;
                nodes[eav->head].supplyWithThetas[0] = eav->mult;
            }
        } else {
            if (eav->tail == eav->head) {
                nodes[eav->tail].supplyWithThetas[0] = 1.0 - eav->mult;
            } else {
                nodes[eav->tail].supplyWithThetas[0] = 1.0;
                nodes[eav->head].supplyWithThetas[0] = -1.0 * eav->mult;
            }
        }
    } else { // var->varType == EQF_VAR
        vector<double> &nodeVals = dynamic_cast<EqFlowVar *>(var)->nodeVals;
        if (var->setLoc == L_var) {
            for (int i = 0; i < numNodes; ++i) {
                nodes[i].supplyWithThetas.resize(s + 1, 0.0);
                nodes[i].supplyWithThetas[0] = -1.0 * nodeVals[i]; 
            } 
        }else {
            for (int i = 0; i < numNodes; ++i) {
                nodes[i].supplyWithThetas.resize(s + 1, 0.0);
                nodes[i].supplyWithThetas[0] = nodeVals[i]; 
            } 
        }
    }
    return;
}

// Running time: O(n)
void GenNetEq::setSupplyWithEqs(Variable *var)
{
    if (var->varType == ARC_VAR) {
        ArcVar *eav = dynamic_cast<ArcVar *>(var);
        // Update the supplies at the given arc
        if (eav->setLoc == L_var) { // eav at LB
            if (eav->tail == eav->head) {
                nodes[eav->tail].supplyWithEq.a0 = -1.0 * (1.0 - eav->mult);
            } else {
                nodes[eav->tail].supplyWithEq.a0 = -1.0;
                nodes[eav->head].supplyWithEq.a0 = eav->mult;
            }
        } else {
            if (eav->tail == eav->head) {
                nodes[eav->tail].supplyWithEq.a0 = 1.0 - eav->mult;
            } else {
                nodes[eav->tail].supplyWithEq.a0 = 1.0;
                nodes[eav->head].supplyWithEq.a0 = -1.0 * eav->mult;
            }
        }
    } else { // var->varType == EQF_VAR
        vector<double> &nodeVals = dynamic_cast<EqFlowVar *>(var)->nodeVals;
        if (var->setLoc == L_var) {
            for (int i = 0; i < numNodes; ++i) {
                nodes[i].supplyWithEq.a0 = -1.0 * nodeVals[i];
            }
        } else { // var->setLoc == U_var
            for (int i = 0; i < numNodes; ++i) {
                nodes[i].supplyWithEq.a0 = nodeVals[i];
            }
        }
    }
    return;
}

/*****************************************************************************/
/* Standard tree computations                                                */
/*****************************************************************************/
// Running time: O(|Ti|)
void GenNetEq::computeFlowsTypeI(TreeTypeI *tI)
{
    ArcVar *extraArc = (*tI).getExtraArc();
    int alpha = extraArc->tail;
    int beta = extraArc->head;
    double muAB = extraArc->mult;
    // x_{alpha,beta} = theta
    extraArc->flowWithEq.a0 = 0.0;
    extraArc->flowWithEq.a1 = 1.0;
    extraArc->flowWithEq.thetaID = THETA;
    // e(alpha) = e(alpha) - theta
    nodes[alpha].supplyWithEq.a1 -= 1.0;  // First term should already be set
    nodes[alpha].supplyWithEq.thetaID = THETA;
    // e(beta) = e(beta) + mu_{alpha,beta} * theta
    nodes[beta].supplyWithEq.a1 += muAB; // First term should already be set
    nodes[beta].supplyWithEq.thetaID = THETA;

    // Perform traversal
    (*tI).resetLeafToRootTraversal();
    int h = (*tI).getRoot();
    int j = (*tI).getPrevThread();
    while (j != h) {
        ArcVar *av = nodes[j].predArc;
        if (av->tail == j) {
            int i = av->head;
            double muJI = av->mult;
            // x_{j,i} = e(j)
            av->flowWithEq.a0 = nodes[j].supplyWithEq.a0;
            av->flowWithEq.a1 = nodes[j].supplyWithEq.a1;
            av->flowWithEq.thetaID = nodes[j].supplyWithEq.thetaID;
            // e(i) += e(j) mu_{j,i}
            nodes[i].supplyWithEq.a0 += nodes[j].supplyWithEq.a0 * muJI;
            nodes[i].supplyWithEq.a1 += nodes[j].supplyWithEq.a1 * muJI;
            nodes[i].supplyWithEq.thetaID = nodes[j].supplyWithEq.thetaID;
        } else { // (av->head == j
            int i = av->tail;
            double muIJ = av->mult;
            // x_{j,i} = -1.0 * e(j) / muIJ
            av->flowWithEq.a0 = nodes[j].supplyWithEq.a0 / (-1.0 * muIJ);
            av->flowWithEq.a1 = nodes[j].supplyWithEq.a1 / (-1.0 * muIJ);
            av->flowWithEq.thetaID = nodes[j].supplyWithEq.thetaID;
            // e(i) += e(j) / mu_{i,j}
            nodes[i].supplyWithEq.a0 += nodes[j].supplyWithEq.a0 / muIJ;
            nodes[i].supplyWithEq.a1 += nodes[j].supplyWithEq.a1 / muIJ;
            nodes[i].supplyWithEq.thetaID = nodes[j].supplyWithEq.thetaID;
        }
        // Reset supplyWithEqs values
        nodes[j].supplyWithEq.a0 = 0.0;
        nodes[j].supplyWithEq.a1 = 0.0;
        nodes[j].supplyWithEq.thetaID = NULL_THETA;
        j = (*tI).getPrevThread();
    }
    // Solve for theta, then reset supplyWithEqs values for root
    double theta = (-1.0 * nodes[h].supplyWithEq.a0)/(nodes[h].supplyWithEq.a1);
    nodes[h].supplyWithEq.a0 = 0.0;
    nodes[h].supplyWithEq.a1 = 0.0;
    nodes[h].supplyWithEq.thetaID = NULL_THETA;

    // NOTE: We always place the new flow value in deltaFlow. If we're 
    // pivoting, this is where we want it. If we're not pivoting, then we 
    // need to copy it over to the flow field, but that's fine since we don't 
    // do this that frequently. 
    vector<ArcVar *> *tIarcs = (*tI).getArcsInTree();
    for (int tAVInd = 0; tAVInd < (*tIarcs).size(); ++tAVInd) {
        ArcVar *av = (*tIarcs)[tAVInd];
        av->deltaFlow = av->flowWithEq.a0 + av->flowWithEq.a1 * theta;
        av->flowWithEq.a0 = 0.0;
        av->flowWithEq.a1 = 0.0;
        av->flowWithEq.thetaID = NULL_THETA;
        checkVarForBlocking(av); // NEW
    }
    // Also need to solve for the flow on the extra arc
    extraArc->deltaFlow = extraArc->flowWithEq.a0 
                        + extraArc->flowWithEq.a1 * theta;
    extraArc->flowWithEq.a0 = 0.0;
    extraArc->flowWithEq.a1 = 0.0;
    extraArc->flowWithEq.thetaID = NULL_THETA;
    checkVarForBlocking(extraArc); // NEW

    return;
}

// Running time: O(s|Ti|)
int GenNetEq::computeFlowsTypeII(TreeTypeII *tII)
{
    int s = setBeqf.size();
    (*tII).resetLeafToRootTraversal();
    int h = (*tII).getRoot();
    int j = (*tII).getPrevThread();
    while (j != h) {
        ArcVar *av = nodes[j].predArc;
        av->flowWithThetas.resize(s + 1, 0.0);
        if (av->tail == j) {
            int i = av->head;
            double muJI = av->mult;
            for (int k = 0; k <= s; ++k) {
                double e_j_k = nodes[j].supplyWithThetas[k]; 
                av->flowWithThetas[k] = e_j_k;
                nodes[i].supplyWithThetas[k] += e_j_k * muJI;
            }
        } else { // (av->head == j
            int i = av->tail;
            double muIJ = av->mult;
            for (int k = 0; k <= s; ++k) {
                double e_j_k = nodes[j].supplyWithThetas[k];
                av->flowWithThetas[k] = (-1.0 * e_j_k) / muIJ;  
                nodes[i].supplyWithThetas[k] += e_j_k / muIJ;
            }
        }
        j = (*tII).getPrevThread();
    }
    return h;
}

