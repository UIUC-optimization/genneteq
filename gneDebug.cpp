/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-07-15                                                          */
/* File: gneDebug.cpp                                                        */
/* Description:                                                              */
/*   Contains implementation details for debugging the GenNetEq algorithm.   */
/*****************************************************************************/
#include "main.h"
#include "matrix.h"
#include "graph.h"
#include "variables.h"
#include "trees.h"
#include "genneteq.h"

#include <gsl/gsl_linalg.h>

#include <cstdio>
#include <cstdlib>

#include <vector>
using std::vector;

// NOTE: Currently deprecated

// Helpful typedef's
typedef vector<ArcVar *>::const_iterator cArcIter;
typedef vector<EqFlowVar *>::const_iterator cEqfIter;

/*****************************************************************************/
/* Debugging functions with pretty-print options                             */
/*****************************************************************************/
void GenNetEq::printAll(bool pretty) const
{
    printf("*** Values computed by GenNetEq ***\n");
    printBasis(pretty);
    printBasicFlow(pretty);
    printPotentials(pretty);
    printReducedCosts(pretty);
    return;
}

void GenNetEq::printBasis(bool pretty) const
{
    printf("Basis is:\n");

    // Print elements of B
    printf("B:=[");
    for (cArcIter xbI = setBarc.begin(); xbI != setBarc.end(); ++xbI) {
        ArcVar *xb = *xbI;
        if (xb->isArtificial) {
            continue;
        }
        printf("(%d,%d),", xb->tail + (pretty ? NODE_IND_OFFSET : 0), 
                           xb->head + (pretty ? NODE_IND_OFFSET : 0));
    }
    for (cEqfIter xbI = setBeqf.begin(); xbI != setBeqf.end(); ++xbI) {
        EqFlowVar *xb = *xbI;
        printf("R_{%d},", xb->eqFlowIndex + (pretty ? EQF_IND_OFFSET : 0));
    }
    printf("]\n");

    // Print elements of L
    printf("L:=[");
    for (cArcIter xlI = setLarc.begin(); xlI != setLarc.end(); ++xlI) {
        ArcVar *xl = *xlI;
        if (xl->isArtificial) {
            continue;
        }
        printf("(%d,%d),", xl->tail + (pretty ? NODE_IND_OFFSET : 0), 
                           xl->head + (pretty ? NODE_IND_OFFSET : 0));
    }
    for (cEqfIter xlI = setLeqf.begin(); xlI != setLeqf.end(); ++xlI) {
        EqFlowVar *xl = *xlI;
        printf("R_{%d},", xl->eqFlowIndex + (pretty ? EQF_IND_OFFSET : 0));
    }
    printf("]\n");

    // Print elements of U
    printf("U:=[");
    for (cArcIter xuI = setUarc.begin(); xuI != setUarc.end(); ++xuI) {
        ArcVar *xu = *xuI;
        if (xu->isArtificial) {
            continue;
        }
        printf("(%d,%d),", xu->tail + (pretty ? NODE_IND_OFFSET : 0), 
                           xu->head + (pretty ? NODE_IND_OFFSET : 0));
    }
    for (cEqfIter xuI = setUeqf.begin(); xuI != setUeqf.end(); ++xuI) {
        EqFlowVar *xu = *xuI;
        printf("R_{%d},", xu->eqFlowIndex + (pretty ? EQF_IND_OFFSET : 0));
    }
    printf("]\n");

    return;
}

void GenNetEq::printFlow(bool pretty) const
{
    printf("Current Flow:\n");
    for (int arcVarInd = 0; arcVarInd < arcVars.size(); ++arcVarInd) {
        ArcVar *av = arcVars[arcVarInd];
        if (av->isArtificial) {
        //if (av->tail == av->head) {
            continue;
        }
        printf("  x_{%d,%d} = %lf\n", 
                av->tail + (pretty ? NODE_IND_OFFSET : 0), 
                av->head + (pretty ? NODE_IND_OFFSET : 0), av->flow);
    }
    for (int eqfVarInd = 0; eqfVarInd < eqFlowVars.size(); ++eqfVarInd) {
        EqFlowVar *eqfv = eqFlowVars[eqfVarInd];
        printf("  x_{%d} = %lf\n", 
                eqfv->eqFlowIndex + (pretty ? EQF_IND_OFFSET : 0), eqfv->flow);
    }
    return;
}

void GenNetEq::printFlow(vector<double> *flowVec, bool pretty) const
{
    printf("y Flow:\n");
    for (int arcVarInd = 0; arcVarInd < arcVars.size(); ++arcVarInd) {
        ArcVar *av = arcVars[arcVarInd];
        if (av->isArtificial) {
        //if (av->tail == av->head) {
            continue;
        }
        printf("  x_{%d,%d} = %lf\n", 
            av->tail + (pretty ? NODE_IND_OFFSET : 0), 
            av->head + (pretty ? NODE_IND_OFFSET : 0), 
            (*flowVec)[av->varIndex]);
    }
    for (int eqfVarInd = 0; eqfVarInd < eqFlowVars.size(); ++eqfVarInd) {
        EqFlowVar *eqfv = eqFlowVars[eqfVarInd];
        printf("  x_{%d} = %lf\n", 
            eqfv->eqFlowIndex + (pretty ? EQF_IND_OFFSET : 0), 
            (*flowVec)[eqfv->varIndex]);
    }
    return;
}

void GenNetEq::printBasicFlow(bool pretty) const
{
    printf("Current Flow for Basic Variables:\n");
    for (cArcIter xbI = setBarc.begin(); xbI != setBarc.end(); ++xbI) {
        ArcVar *xb = *xbI;
        printf("  x_{%d,%d}", xb->tail + (pretty ? NODE_IND_OFFSET : 0), 
                              xb->head + (pretty ? NODE_IND_OFFSET : 0));
        if (xb->isArtificial) {
            printf("-F");
        }
        printf(" = %lf\n", xb->flow);
    }
    for (cEqfIter xbI = setBeqf.begin(); xbI != setBeqf.end(); ++xbI) {
        EqFlowVar *xb = *xbI;
        printf("  x_{%d} = %lf\n", 
            xb->eqFlowIndex + (pretty ? EQF_IND_OFFSET : 0), 
            xb->flow);
    }
    return;
}

void GenNetEq::printPotentials(bool pretty) const
{
    printf("Current Potentials:\n");
    for (int i = 0; i < numNodes; ++i) {
        printf("  pi_{%d} = %lf\n", 
            i + (pretty ? NODE_IND_OFFSET : 0), nodes[i].potential);
    }
    return;
}

void GenNetEq::printReducedCosts(bool pretty) const
{
    printf("Current Reduced Costs:\n");
    for (int arcVarInd = 0; arcVarInd < arcVars.size(); ++arcVarInd) {
        ArcVar *av = arcVars[arcVarInd];
        if (av->isArtificial) {
        //if (av->tail == av->head) {
            continue;
        }
        printf("  c^pi_{%d,%d} = %lf\n", 
            av->tail + (pretty ? NODE_IND_OFFSET : 0), 
            av->head + (pretty ? NODE_IND_OFFSET : 0), av->redCost);
    }
    for (int eqfVarInd = 0; eqfVarInd < eqFlowVars.size(); ++eqfVarInd) {
        EqFlowVar *eqfv = eqFlowVars[eqfVarInd];
        printf("  c^pi_{%d} = %lf\n", 
            eqfv->eqFlowIndex + (pretty ? EQF_IND_OFFSET : 0), eqfv->redCost);
    }
    return;
}

void GenNetEq::printSupply(bool pretty) const
{
    printf("Current Supply:\n");
    for (int i = 0; i < numNodes; ++i) {
        printf("  b(%d) = %lf\n", 
            i + (pretty ? NODE_IND_OFFSET : 0), (*supply)[i]);
    }
    return;
}

/*****************************************************************************/
/* Miscellaneous debugging functions                                         */
/*****************************************************************************/
void GenNetEq::printFlowsWithThetas() const
{
    printf("Flows with variables:\n");
    for (int arcVarInd = 0; arcVarInd < arcVars.size(); ++arcVarInd) {
        ArcVar *av = arcVars[arcVarInd];
        printf("  x_{%d,%d}", av->tail, av->head);
        if (av->isArtificial) {
            printf("-F");
        }
        printf(" = "); 
        printInTermsOfThetas(&(av->flowWithThetas));
    }
    for (int eqfVarInd = 0; eqfVarInd < eqFlowVars.size(); ++eqfVarInd) {
        EqFlowVar *eqfv = eqFlowVars[eqfVarInd];
        printf("  x_{%d} = ", eqfv->eqFlowIndex);
        printInTermsOfThetas(&(eqfv->flowWithThetas));
    }
    return;
}

void GenNetEq::printSupplyWithThetas() const
{
    printf("Supply with variables:\n");
    for (int i = 0; i < numNodes; ++i) {
        printf("  b^(%d) = ", i);
        printInTermsOfThetas(&(nodes[i].supplyWithThetas));
    }
    return;
}

inline
void GenNetEq::printInTermsOfThetas(vector<double> *thetaVals) const
{
    printf("[");
    if (thetaVals == NULL) {
        printf("NULL]\n");
    } else {
        for (int k = 0; k < (*thetaVals).size(); ++k) {
            printf("%lf,", (*thetaVals)[k]);
        }
        printf("]\n");
    }
    fflush(stdout);
    return;
}

void GenNetEq::printFlowsWithEqs() const
{
    printf("Flows with equations:\n");
    for (int arcVarInd = 0; arcVarInd < arcVars.size(); ++arcVarInd) {
        ArcVar *av = arcVars[arcVarInd];
        printf("  x_{%d,%d}", av->tail, av->head);
        if (av->isArtificial) {
            printf("-F");
        }
        printf(" = ");
        printInTermsOfThetaEquation(&(av->flowWithEq));
    }
    for (int eqfVarInd = 0; eqfVarInd < eqFlowVars.size(); ++eqfVarInd) {
        EqFlowVar *eqfv = eqFlowVars[eqfVarInd];
        printf("  x_{%d} = ", eqfv->eqFlowIndex);
        printInTermsOfThetaEquation(&(eqfv->flowWithEq));
    }
    return;
}

void GenNetEq::printSupplyWithEqs() const
{
    printf("Supply with equations:\n");
    for (int i = 0; i < numNodes; ++i) {
        printf("  b^(%d) = ", i);
        printInTermsOfThetaEquation(&(nodes[i].supplyWithEq));
    }
    return;
}

void GenNetEq::printPotentialsWithEqs() const
{
    printf("Potentials with equations:\n");
    for (int i = 0; i < numNodes; ++i) {
        printf("  pi_{%d} = ", i);
        printInTermsOfThetaEquation(&(nodes[i].potentialWithEq));
    }
    return;
}

inline
void GenNetEq::printInTermsOfThetaEquation(ThetaEquation *thetaEq) const
{
    printf("%lf + %lf Theta_{%d}\n", thetaEq->a0, thetaEq->a1, thetaEq->thetaID);

//    if (thetaEq->thetaID == NULL_THETA) {
//        printf("\n");
//    } else {
//        printf(" + %lf Theta", thetaEq->a1);
//        if (thetaEq->thetaID == THETA) {
//            printf("\n");
//        } else {
//            printf("_{%d}\n", thetaEq->thetaID);
//        }
//    }
    return;
}

void GenNetEq::printGSLmatrix(gsl_matrix *mat, int m, int n) const
{
    printf("Matrix := \n");
    for (int i = 0; i < m; ++i) {
        printf("  |");
        for (int j = 0; j < n; ++j) {
            printf(" %lf", gsl_matrix_get(mat, i, j)); 
        }
        printf(" |\n");
    }
    return;
}

void GenNetEq::printGSLvector(gsl_vector *x, int n) const
{
    printf("Vector := {");
    for (int i = 0; i < n; ++i) {
        printf("%lf,", gsl_vector_get(x, i));
    }
    printf("}\n");
    return;
}

/*****************************************************************************/
/* Debugging functions for exporting basis in Octave format                  */
/*****************************************************************************/
void GenNetEq::printComps(int iterNum, bool pretty) const
{
    char octaveFilename[33];
    char genneteqFilename[32];

    sprintf(octaveFilename, "../debug/octaveRaw/iter%d.m", iterNum);
    sprintf(genneteqFilename, "../debug/genneteq/iter%d.g", iterNum);

    exportBasisInOctave(octaveFilename, pretty);
    exportBasisInGenNetEq(genneteqFilename, pretty);

    return;
}

void GenNetEq::exportBasisInOctave(const char *filename, bool pretty) const
{
    FILE *outFile;
    outFile = fopen(filename, "w");
    if (outFile == NULL) {
        printf("Unable to open file for Octave export: %s\n", filename);
        exit(-1);
    }

    // Print header information
    fprintf(outFile, "#!/usr/bin/octave\n");

    // Print out the constraint matrix, costs, capacities, and supply vector
    fprintf(outFile, "### Input Parameters\n");
    for (int arcVarInd = 0; arcVarInd < arcVars.size(); ++arcVarInd) {
        ArcVar *av = arcVars[arcVarInd];
        int i = av->tail + (pretty ? NODE_IND_OFFSET : 0);
        int j = av->head + (pretty ? NODE_IND_OFFSET : 0);
        fprintf(outFile, "A_%d_%d = [", i, j);
        for (int k = 0; k < numNodes; ++k) {
            if (k == i) {
                if (i == j) {
                    fprintf(outFile, "%lf;", 1.0 - (av->mult));
                } else {
                    fprintf(outFile, "%lf;", 1.0);
                }
            } else if (k == j) {
                fprintf(outFile, "%lf;", -1.0 * (av->mult));
            } else { 
                fprintf(outFile, "%lf;", 0.0);
            }
        }
        fprintf(outFile, "];\n");
        fprintf(outFile, "c_%d_%d = %lf;\n", i, j, av->cost);
        if (i == j) {
            // Done to prevent unnecessary printing of extra characters
            fprintf(outFile, "u_%d_%d = %lf;\n\n", i, j, 100000000000000.0);
        } else {
            fprintf(outFile, "u_%d_%d = %lf;\n\n", i, j, av->UB);
        }
    }
    for (int eqfVarInd = 0; eqfVarInd < eqFlowVars.size(); ++eqfVarInd) {
        EqFlowVar *eqfv = eqFlowVars[eqfVarInd];
        int eqfIndex = eqfv->eqFlowIndex + (pretty ? EQF_IND_OFFSET : 0);
        fprintf(outFile, "A_%d = [", eqfIndex);
        for (int k = 0; k < numNodes; ++k) {
            //fprintf(outFile, "%lf;", (*(eqfv->nodeVals))[k]); // DEPRECATED
            fprintf(outFile, "%lf;", eqfv->nodeVals[k]);
        }
        fprintf(outFile, "];\n");
        fprintf(outFile, "c_%d = %lf;\n", eqfIndex, eqfv->cost);
        fprintf(outFile, "u_%d = %lf;\n\n", eqfIndex, eqfv->UB);
    }
    fprintf(outFile, "\nb = [");
    for (int i = 0; i < numNodes; ++i) {
        fprintf(outFile, "%lf;", (*supply)[i]);
    }
    fprintf(outFile, "];\n\n");

    // Print out updated supply vector, b^hat, based on upper bounded vars
    fprintf(outFile, "### Updating supply vector based on upper bounds\n");
    fprintf(outFile, "bHat = b");
    for (int nbArcInd = 0; nbArcInd < setNBarc.size(); ++nbArcInd) {
        ArcVar *xu = setNBarc[nbArcInd];
        if (xu->setLoc == L_var) { continue; }
        int i = xu->tail + (pretty ? NODE_IND_OFFSET : 0);
        int j = xu->head + (pretty ? NODE_IND_OFFSET : 0);
        fprintf(outFile, " - (u_%d_%d * A_%d_%d)", i, j, i, j);
    }
    for (int nbEqfInd = 0; nbEqfInd < setNBeqf.size(); ++nbEqfInd) {
        EqFlowVar *xu = setNBeqf[nbEqfInd];
        if (xu->setLoc == L_var) { continue; }
        int eqfIndex = xu->eqFlowIndex + (pretty ? EQF_IND_OFFSET : 0);
        fprintf(outFile, " - (u_%d * A_%d)", eqfIndex, eqfIndex);
    }
    fprintf(outFile, ";\n\n");

    // Print out basis structure
    fprintf(outFile, "### Generating basis matrix\n");
    fprintf(outFile, "B = [");
    for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
        ArcVar *xb = setBarc[bArcInd];
        fprintf(outFile, "A_%d_%d,", xb->tail + (pretty ? NODE_IND_OFFSET :0), 
                                     xb->head + (pretty ? NODE_IND_OFFSET :0));
    }
    for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
        EqFlowVar *xb = setBeqf[bEqfInd];
        fprintf(outFile, "A_%d,", 
            xb->eqFlowIndex + (pretty ? EQF_IND_OFFSET : 0));
    }
    fprintf(outFile, "];\n\n");

    // Print out basis cost vector
    fprintf(outFile, "### Generating basis cost vector\n");
    fprintf(outFile, "c_BT = [");
    for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
        ArcVar *xb = setBarc[bArcInd];
        fprintf(outFile, "c_%d_%d,", xb->tail + (pretty ? NODE_IND_OFFSET :0), 
                                     xb->head + (pretty ? NODE_IND_OFFSET :0));
    }
    for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
        EqFlowVar *xb = setBeqf[bEqfInd];
        fprintf(outFile, "c_%d,", 
            xb->eqFlowIndex + (pretty ? EQF_IND_OFFSET : 0));
    }
    fprintf(outFile, "];\n\n");

    // Compute flow vector for the given basis
    fprintf(outFile, "### Computing flow vector\n");
    fprintf(outFile, "x_B = B \\ bHat\n\n");

    // Compute node potentials for the given basis
    fprintf(outFile, "### Computing node potentials\n");
    fprintf(outFile, "Binv = inv(B);\n");
    fprintf(outFile, "yT = c_BT * Binv;\n");
    fprintf(outFile, "pi = transpose(yT)\n\n");

    // Print out equations to solve for the reduced costs of all variables
    fprintf(outFile, "### Computing reduced costs\n");
    for (int arcVarInd = 0; arcVarInd < arcVars.size(); ++arcVarInd) {
        ArcVar *av = arcVars[arcVarInd];
        int i = av->tail + (pretty ? NODE_IND_OFFSET : 0);
        int j = av->head + (pretty ? NODE_IND_OFFSET : 0);
        fprintf(outFile, "cPi_%d_%d = c_%d_%d - yT * A_%d_%d\n", 
                            i, j, i, j, i, j);
    }
    for (int eqfVarInd = 0; eqfVarInd < eqFlowVars.size(); ++eqfVarInd) {
        EqFlowVar *eqfv = eqFlowVars[eqfVarInd];
        int eqfIndex = eqfv->eqFlowIndex + (pretty ? EQF_IND_OFFSET : 0);
        fprintf(outFile, "cPi_%d = c_%d - yT * A_%d\n", 
                            eqfIndex, eqfIndex, eqfIndex);
    }
    fprintf(outFile, "\n");

    fclose(outFile);
    return;
}

void GenNetEq::exportBasisInGenNetEq(const char *filename, bool pretty) const
{
    FILE *outFile;
    outFile = fopen(filename, "w");
    if (outFile == NULL) {
        printf("Unable to open file for GenNetEq export\n");
        exit(-1);
    }
    // Printing basic flow values
    fprintf(outFile, "x_B = \n\n");
    for (cArcIter xbI = setBarc.begin(); xbI != setBarc.end(); ++xbI) {
        ArcVar *xb = *xbI;
        fprintf(outFile, "  %lf\n", xb->flow);
    }
    for (cEqfIter xbI = setBeqf.begin(); xbI != setBeqf.end(); ++xbI) {
        EqFlowVar *xb = *xbI;
        fprintf(outFile, "  %lf\n", xb->flow);
    }
    // Printing node potentials
    fprintf(outFile, "\npi = \n\n");
    for (int i = 0; i < numNodes; ++i) {
        fprintf(outFile, "  %lf\n", nodes[i].potential);
    }
    fprintf(outFile, "\n"); 
    // Printing reduced costs
    for (int arcVarInd = 0; arcVarInd < arcVars.size(); ++arcVarInd) {
        ArcVar *av = arcVars[arcVarInd];
        int i = av->tail + (pretty ? NODE_IND_OFFSET : 0);
        int j = av->head + (pretty ? NODE_IND_OFFSET : 0);
        fprintf(outFile, "cPi_%d_%d = %lf\n", i, j, av->redCost);
    }
    for (int eqfVarInd = 0; eqfVarInd < eqFlowVars.size(); ++eqfVarInd) {
        EqFlowVar *eqfv = eqFlowVars[eqfVarInd];
        int eqfIndex = eqfv->eqFlowIndex + (pretty ? EQF_IND_OFFSET : 0);
        fprintf(outFile, "cPi_%d = %lf\n", eqfIndex, eqfv->redCost);
    }
    fprintf(outFile, "\n");

    fclose(outFile);
    return;  
}

