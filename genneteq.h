/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-06-16                                                          */
/* File: genneteq.h                                                          */
/* Description:                                                              */
/*   Contains design details for the GenNetEq algorithm.                     */
/*****************************************************************************/
#ifndef GENNETEQ_H
#define GENNETEQ_H

// Required include's
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vector>
using std::vector;
#include <list>
using std::list;
#include <queue>
using std::priority_queue;

#include <gsl/gsl_linalg.h>

// Forward declarations
class Graph;
class Node;
class Variable;
class ArcVar;
class EqFlowVar;
class Tree;
struct PivotVars;
struct VarWithViol;

class GenNetEq
{
  public:
    GenNetEq();
    ~GenNetEq();

    int solve(bool usePresolve);
    void initialize(Graph *g, int initType);

  protected:
    // Variables
    int numNodes;
    int numArcs;
    int numEqualFlowSets;
    double bigM;
    bool presolve;
    int initialType;

    double flowCost;

// Stat Counters
    bool isOptimal;
    bool isInfeasible;
    bool isUnbounded;

    int phaseIiters;
    int phaseIIiters;
    int loopIter;
    int infeasLoopIter;

    int consecutiveCycles;
    int consecutiveDegenPivots;
    int totalCyclesDetected;
    int totalDegenPivots;

// Parameters
    bool trackLowestVarIndex;
    bool selectRandomBlockingVar;
    bool selectRandomOffendingVar;
    bool useRandomJumpWhenCycling;
    int numLargestViolVarsToTrack;
    int majorIterationLength;
    double violationThreshold;

// Data structures and other variables for algorithm
    vector<Node> nodes;

    vector<Variable *> vars;
    vector<ArcVar *> arcVars;
    vector<EqFlowVar *> eqFlowVars;

    vector<ArcVar *> setBarc;
    vector<EqFlowVar *> setBeqf;

    vector<ArcVar *> setNBarc;
    vector<EqFlowVar *> setNBeqf;

    vector<TreeTypeI *> setBtI;
    vector<TreeTypeII *> setBtII;

    // DEPRECATED
//    vector<ArcVar *> setLarc;
//    vector<EqFlowVar *> setLeqf;
//    vector<ArcVar *> setUarc;
//    vector<EqFlowVar *> setUeqf;

    // Maybe move this into the variables, depending on if we need it
    vector<int> useArtificialCosts;

    // Used for pivoting
    double delta;
    Variable *eVar;
    Variable *lVar;
    vector<Variable *> blockingVars;
    vector<Variable *> offendingVars;

    bool cyclingDetected;
    list<PivotVars *> cycleList;
    int lastEVarInd;
	double lastViolation;
	double baselineViolation;

    // Candidate List Variables
    priority_queue<VarWithViol*, vector<VarWithViol*>, violcomp> candidateList;
    int itersSinceUpdate;
    vector<Variable *> candidateList2;

    int clMaxSize;
    int nbArcPivotInd;
    int nbEqfPivotInd;
    int numArcsToAlwaysCheck;
    int numEqfsToAlwaysCheck;

    // Used for initialization procedure
    vector<ArcVar*> initCands;
    vector<ArcVar*> initTrees;
    vector<double> supplyToRoot;
    vector<int> initPreds;
    vector<int> initComps;
    int currentCompNum;
    int currNumInfeas;

    // Used for stats tracking
	int totalPivots;
	int arcArcPivots;
	int arcEFSPivots;
	int EFSArcPivots;
	int EFSEFSPivots;
	int sameTIPivots;
	int sameTIIPivots;
	int tItIPivots;
	int tIItIIPivots;
	int tItIIPivots;
	int EFSPivots;

	float avgTIsize;
	float avgTIIsize;

    // Functions
// Miscellaneous
    void initializeStats(bool usePresolve);
    void printTerminationStatus();
    void checkFeasibility();

    void computeFlowCost();
    void setVariableCosts();

    double computeArcRedCost(ArcVar *av);                           // inline
    double computeEqfRedCost(EqFlowVar *eqfv);                      // inline
    void checkVarForBlocking(Variable *var);                        // inline
    double computeDeltaVar(Variable *var);                          // inline
    bool checkVarFlow(Variable *var, bool useDeltaFlow = false);    // inline

// Init
    void initializeNodes(Graph *g);
    void initializeVars(Graph *g);
    void formInitialBFS(Graph *g);

    void addInitialArc(ArcVar *arc);
    int computeNumInfeasible();
    double computeNewInfeasible(ArcVar* arc);
    void computeInitFlowOnArc(int curr, int pred, double delta);
    void computeInitFlowInfeasible(int curr, int pred, double delta, 
            bool useFlow, vector<double>& strCopy, vector<int>& infeasArcs);
    double score(ArcVar* arc, double su, double sv, bool forward);

    ArcVar* findAV(int u, int v);
    ArcVar* findBasicAV(int u, int v);

// Flow
    void computeFlows(int updateOption = 0);
    void setFlowValues(int updateOption);
    void computeFlowsPivotingTI(ArcVar *eav, TreeTypeI *tU, TreeTypeI *tV);
    void computeFlowsPivotingTII(Variable *var);

    // Used with both flow computations
    void computeFlowsWithEqFlowSets();
    void computeNumericFlowsAndSupply(gsl_vector *x);

    // Used with standard flow computation 
    void computeBoundedFlows(bool useEqFlowSets);

    void setSupplyWithThetas();
    void setSupplyWithEqs();
    void setSupplyWithThetas(Variable *var);
    void setSupplyWithEqs(Variable *var);

    // Used to compute the flows on the type I and type II trees
    void computeFlowsTypeI(TreeTypeI *tI);
    int computeFlowsTypeII(TreeTypeII *tII);

// Potential
    void computePotentials();
    void computePotentialsInTermsOfRoot(Tree *tX);
    void finalizePotentialsTypeI(TreeTypeI *tI);
    void computePotentialsWithEqFlowSets();
    void computeReducedCosts(bool includeBasicVars = false);
    void verifyOptimality();

// Pivot
    void pivot();
    void pivotTI(ArcVar *eav);
    void pivotTII();

    void identifyLeavingVar();

    void updateArtificialCosts();
    void updateBasis();
    void checkForCycling();
    void updatePivotingStats();

    double computeMinDelta(Variable *var);
    void convertBasicArcsToTrees();

    void insertIntoB(Variable *var);
    void removeFromB(Variable *var);
    void insertIntoL(Variable *var);
    void removeFromL(Variable *var);
    void insertIntoU(Variable *var);
    void removeFromU(Variable *var);

// Pivot Rules
    void identifyEnteringVariable();
    Variable *selectEnteringVariable();

    Variable *getVarWithLargestViolation(bool useWeighted = false);
    Variable *getVarWithLargestViolationArcsFirst();
    Variable *getVarWithLargestViolationEqfsFirst();

    Variable *getFirstVarWithViolation();
    Variable *getVarWithSteepestDescent();

    Variable *getRandomVarWithViolation();
    Variable *getRandomVarWeightedByViolation();
    Variable *getRandomVarWithLargeViolation(bool useWeighted = false);

	Variable* getVarBySimAnnealingArcsFirst();
	Variable* getVarBySimAnnealingEqfsFirst();

    // Helper Functions
    void computeArcsWithLargestViolation(double *largestViolation);
    void computeEqfsWithLargestViolation(double *largestViolation, 
                                         bool useWeighted = false);
    void checkVarForOffending(Variable *var, double varViolation,   // inline
                              double *largestViolation);

    // Candidate List Pivoting Functions
    Variable *getNextCandidate(bool useWeighted = false);
    void updateCandidateList(bool useWeighted);
    void clearCandidateList();

    Variable *getNextCandidate2(bool useWeighted = false);
    void updateCandidateList2();
    void tryToAddArc();
    void tryToAddEqf();

    Variable *quickSelect(bool useWeighted = false);

public:
/*
// Debug
    void printAll(bool pretty = false) const;
    void printBasis(bool pretty = false) const;
    void printFlow(bool pretty = false) const;
    void printFlow(vector<double> *flowVec, bool pretty = false) const;
    void printBasicFlow(bool pretty = false) const;
    void printPotentials(bool pretty = false) const;
    void printReducedCosts(bool pretty = false) const;
    void printSupply(bool pretty = false) const;

    // Printing contents of workspace data structures
    void printFlowsWithThetas() const;
    void printSupplyWithThetas() const;
    void printInTermsOfThetas(vector<double> *thetaVals) const;

    void printFlowsWithEqs() const;
    void printSupplyWithEqs() const;
    void printPotentialsWithEqs() const;
    void printInTermsOfThetaEquation(ThetaEquation *thetaEq) const;

    void printGSLmatrix(gsl_matrix *mat, int m, int n) const;
    void printGSLvector(gsl_vector *x, int n) const;

    // Exporting basis to Octave for comparison
    void printComps(int iterNum, bool pretty = false) const;
    void exportBasisInOctave(const char *filename, bool pretty = false) const;
    void exportBasisInGenNetEq(const char *filename, bool pretty=false) const;
// */

  private:
    // Nothing
};

/*****************************************************************************/
/* Inline Function Definitions                                               */
/*****************************************************************************/
// Running time: O(1)
inline
double GenNetEq::computeArcRedCost(ArcVar *av)
{
    av->redCost = av->cost - nodes[av->tail].potential 
                + av->mult * (nodes[av->head].potential);
    return av->redCost;
}

// Running time: O(n)
inline
double GenNetEq::computeEqfRedCost(EqFlowVar *eqfv) 
{
    double eqfvRC = eqfv->cost;
    for (int i = 0; i < numNodes; ++i) {
        eqfvRC -= ((eqfv->nodeVals[i]) * (nodes[i].potential));
    }
    eqfv->redCost = eqfvRC;
    return eqfvRC;
}

/*****************************************************************************/
// Running time: O(1)
inline 
void GenNetEq::checkVarForBlocking(Variable *var)
{
    double deltaVar = computeDeltaVar(var);
    // DEBUG:
//    if (deltaVar < 0) {
//        printf("ERROR (%d): delta for variable %d is negative: %.16lf\n", 
//            loopIter, var->varIndex, deltaVar);
//        printf("Var %d: xFlow: %.16lf: yFlow: %.16lf\n", var->varIndex, 
//            var->flow, var->deltaFlow);
//    }

    if (delta > deltaVar + ROUNDOFF_TOLERANCE) {
        blockingVars.clear();
    }
    if (delta >= deltaVar - ROUNDOFF_TOLERANCE) {
        delta = deltaVar;
//        if ((!trackLowestVarIndex) || (blockingVars.empty()) ||
//            (blockingVars.front()->varIndex < var->varIndex)) {
            blockingVars.push_back(var);
//        } else {
//            blockingVars.push_back(blockingVars.front());
//            blockingVars[0] = var;
//        }
    }
    return;
} 

// Running time: O(1)
inline 
double GenNetEq::computeDeltaVar(Variable *var)
{
    // Compute the maximum amount of flow that can be pushed over the entering
    // edge before the given basic variable reaches its upper or lower bound.
    double xVar = var->flow;
    double yVar = var->deltaFlow;

    if (yVar > ROUNDOFF_TOLERANCE) {
        // If the flow on this variable increases, this stops when it hits its 
        // upper bound, unless it's already above it; in that case, it can 
        // increase as much as it wants.
        if ((initialType == INFEASIBLE) && presolve && 
            (useArtificialCosts[var->varIndex] > 0)) {
            return std::numeric_limits<double>::max();
        } else {
            if (fabs(var->UB - xVar) < ROUNDOFF_TOLERANCE) {
                return 0.0;
            } else {
                return (var->UB - xVar) / yVar;
            }
        }
    } else if (yVar < (-1.0 * ROUNDOFF_TOLERANCE)) {
        // Otherwise, the flow on the edge decreases, so we stop when we hit 
        // the lower bound (again, unless the flow is negative, in which case 
        // it can go as far down as it likes).
        if ((initialType == INFEASIBLE) && presolve && 
            (useArtificialCosts[var->varIndex] < 0)) {
            return std::numeric_limits<double>::max();
        } else {
            if (fabs(xVar) < ROUNDOFF_TOLERANCE) {
                return 0.0;
            } else {
                return xVar / (-1.0 * yVar);
            }
        }
    } else {
        // Otherwise, variable doesn't change, so it won't block
        return std::numeric_limits<double>::max();
    }
}

/*****************************************************************************/
// Running time: O(1)
inline
bool GenNetEq::checkVarFlow(Variable *var, bool useDeltaFlow)
{
    double varFlow = useDeltaFlow ? var->deltaFlow : var->flow;
    if (((initialType != INFEASIBLE) || 
         (useArtificialCosts[var->varIndex] == 0.0)) &&
        ((varFlow < var->LB - ROUNDOFF_TOLERANCE) ||
         (varFlow > var->UB + ROUNDOFF_TOLERANCE))) {
        printf("ERROR (%d): Var %d's ", loopIter, var->varIndex); 
        if (useDeltaFlow) { printf("recomputed"); }
        else { printf("current"); }
        printf(" flow is not within bounds: %.16lf != [%.16lf,%.16lf]\n",  
                varFlow, var->LB, var->UB);
        return false;
    }
    return true;
}

#endif // GENNETEQ_H

