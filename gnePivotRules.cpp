/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-12-20                                                          */
/* File: gnePivotRules.cpp                                                   */
/* Description:                                                              */
/*   Contains implementation details for the pivoting rules used by the      */
/*   GenNetEq algorithm.                                                     */
/*****************************************************************************/
#include "main.h"
#include "variables.h"
#include "trees.h"
#include "genneteq.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vector>
using std::vector;

extern int phaseIpivot;
extern int phaseIIpivot;

/*****************************************************************************/
/* Pivoting Rules Protected Functions                                        */
/*****************************************************************************/
// Running time: O(t_select)
void GenNetEq::identifyEnteringVariable()
{
    eVar = selectEnteringVariable();
    //eVar = getNextCandidate2(true);
    return;
}

// Running time: O(t_select)
Variable *GenNetEq::selectEnteringVariable()
{
    // NOTE: For running the tests, I think we should modify the code to just 
    // do the following instead of executing the switch statement:
    // return <appropriatePivotRule>(); 

    // If cycling, return random variable with large (weighted) violation.
    if (cyclingDetected && useRandomJumpWhenCycling) {
//        printf("Cycle detected, applying random pivot\n");
        return getRandomVarWithLargeViolation();
    } // else apply standard pivoting rules

    if (presolve) { // Pivoting rule to use in Phase I
        switch (phaseIpivot) {
          case LV:
            return getVarWithLargestViolation();
          case LVW:
            return getVarWithLargestViolation(true);
          case LVA:
            return getVarWithLargestViolationArcsFirst();
		  case LVE:
			return getVarWithLargestViolationEqfsFirst();
          case FV:
            return getFirstVarWithViolation();
          case SD:
            return getVarWithSteepestDescent();
          case RV:
            return getRandomVarWithViolation();
          case RWV:
            return getRandomVarWeightedByViolation();
          case RLV:
            return getRandomVarWithLargeViolation();
          case RLVW:
            return getRandomVarWithLargeViolation(true);
          case CL:
            return getNextCandidate();
          case CLW:
            return getNextCandidate(true);
		  case SAA:
			return getVarBySimAnnealingArcsFirst();
		  case SAE:
			return getVarBySimAnnealingEqfsFirst();
          case CL2:
            return getNextCandidate2();
          case CLW2:
            return getNextCandidate2(true);
          case QS:
            return quickSelect();
          case QSW:
            return quickSelect(true);
        }
    } else { // Pivoting rule to use in Phase II
        switch (phaseIIpivot) {
          case LV:
            return getVarWithLargestViolation();
          case LVW:
            return getVarWithLargestViolation(true);
          case LVA:
            return getVarWithLargestViolationArcsFirst();
		  case LVE:
			return getVarWithLargestViolationEqfsFirst();
          case FV:
            return getFirstVarWithViolation();
          case SD:
            return getVarWithSteepestDescent();
          case RV:
            return getRandomVarWithViolation();
          case RWV:
            return getRandomVarWeightedByViolation();
          case RLV:
            return getRandomVarWithLargeViolation();
          case RLVW:
            return getRandomVarWithLargeViolation(true);
          case CL:
            return getNextCandidate();
          case CLW:
            return getNextCandidate(true);
		  case SAA:
			return getVarBySimAnnealingArcsFirst();
		  case SAE:
			return getVarBySimAnnealingEqfsFirst();
          case CL2:
            return getNextCandidate2();
          case CLW2:
            return getNextCandidate2(true);
          case QS:
            return quickSelect();
          case QSW:
            return quickSelect(true);
        }
    }
}

/*****************************************************************************/
/* Pivoting Rules                                                            */
/*****************************************************************************/
// Running time: O(m' + np)
Variable *GenNetEq::getVarWithLargestViolation(bool useWeighted)
{
    isOptimal = true;
    offendingVars.clear();
    double largestViolation = -1.0;

    computeArcsWithLargestViolation(&largestViolation);
    computeEqfsWithLargestViolation(&largestViolation, useWeighted);
	lastViolation = largestViolation;

    // DEBUG:
//    printf("---Offendinging (%lf): ", largestViolation);
//    for (int i = 0; i < offendingVars.size(); ++i) {
//        printf("%d,", offendingVars[i]->varIndex);
//    }
//    printf("\n");

    if (offendingVars.empty()) {
        return NULL;
    } else if (selectRandomOffendingVar) {
        return offendingVars[drand48() * offendingVars.size()];
    } else {
        return offendingVars.front();
    }
}

// Running time: O(m' + np)
Variable *GenNetEq::getVarWithLargestViolationArcsFirst()
{
    isOptimal = true;
    offendingVars.clear();
    double largestViolation = -1.0;

    computeArcsWithLargestViolation(&largestViolation);
	lastViolation = largestViolation;
    if (offendingVars.empty()) {
        computeEqfsWithLargestViolation(&largestViolation, false);
		lastViolation = largestViolation;
    }

    if (offendingVars.empty()) {
        return NULL;
    } else if (selectRandomOffendingVar) {
        return offendingVars[drand48() * offendingVars.size()];
    } else {
        return offendingVars.front();
    }
}

// Running time: O(m' + np)
Variable *GenNetEq::getVarWithLargestViolationEqfsFirst()
{
    isOptimal = true;
    offendingVars.clear();
    double largestViolation = -1.0;

    computeEqfsWithLargestViolation(&largestViolation, false);
	lastViolation = largestViolation;
    if (offendingVars.empty()) {
        computeArcsWithLargestViolation(&largestViolation);
		lastViolation = largestViolation;
    }

    if (offendingVars.empty()) {
        return NULL;
    } else if (selectRandomOffendingVar) {
        return offendingVars[drand48() * offendingVars.size()];
    } else {
        return offendingVars.front();
    }
}

// Running time: O(m' + np)
Variable *GenNetEq::getFirstVarWithViolation()
{
    isOptimal = true;
    for (int nbArcInd = 0; nbArcInd < setNBarc.size(); ++nbArcInd) {
        ArcVar *nba = setNBarc[nbArcInd];
        double rc = computeArcRedCost(nba);
        if (((nba->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((nba->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            return nba;
        }
    }
    for (int nbEqfInd = 0; nbEqfInd < setNBeqf.size(); ++nbEqfInd) {
        EqFlowVar *nbeq = setNBeqf[nbEqfInd];
        double rc = computeEqfRedCost(nbeq);
        if (((nbeq->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((nbeq->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            return nbeq;
        }
    }
    return NULL;
}

// Running time: O((m' + p)(ns + s^2 + g(s)))
Variable *GenNetEq::getVarWithSteepestDescent()
{
    isOptimal = true;
    offendingVars.clear();
    double bestImprovement = -1.0;

    // Now go through all non-basic variables and compute their reduced costs. 
    // If a variable violates its optimality conditions, then compute the 
    // minimum delta value to see what the total change in cost will be if 
    // pivoting on this variable. [O((m' + p)(ns + s^2 + g(s)))]
    for (int nbArcInd = 0; nbArcInd < setNBarc.size(); ++nbArcInd) {
        ArcVar *nba = setNBarc[nbArcInd];
        double rc = computeArcRedCost(nba);
        if (((nba->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((nba->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            double varDelta = computeMinDelta(nba);
            double varImprovement = varDelta * fabs(rc);
            checkVarForOffending(nba, varImprovement, &bestImprovement);
        }
    }
    for (int nbEqfInd = 0; nbEqfInd < setNBeqf.size(); ++nbEqfInd) {
        EqFlowVar *nbeq = setNBeqf[nbEqfInd];
        double rc = computeEqfRedCost(nbeq);
        if (((nbeq->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((nbeq->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            double varDelta = computeMinDelta(nbeq);
            double varImprovement = varDelta * fabs(rc);
            checkVarForOffending(nbeq, varImprovement, &bestImprovement);
        }
    }
    if (offendingVars.empty()) {
        return NULL;
    } else if (selectRandomOffendingVar) {
        return offendingVars[drand48() * offendingVars.size()];
    } else {
        return offendingVars.front();
    }
}

// Running time: O(m' + np)
Variable *GenNetEq::getRandomVarWithViolation()
{
    computeReducedCosts();
    if (offendingVars.empty()) {
        return NULL;
    } else {
        return offendingVars[drand48() * offendingVars.size()];
    }
}

// Running time: O(m' + np)
Variable *GenNetEq::getRandomVarWeightedByViolation()
{
    computeReducedCosts();
    if (!offendingVars.empty()) {
        double violationSum = 0.0;
        for (int i = 0; i < offendingVars.size(); ++i) {
            Variable *var = offendingVars[i];
            violationSum += fabs(var->redCost);
        }
        double violationLevel = drand48() * violationSum;
        double runningViolationSum = 0.0;
        for (int i = 0; i < offendingVars.size(); ++i) {
            Variable *var = offendingVars[i];
            runningViolationSum += fabs(var->redCost);
            if (violationLevel <= runningViolationSum) {
                return var;
            }
        }
    }
    return NULL;
}

// Running time: O(m' + np)
Variable *GenNetEq::getRandomVarWithLargeViolation(bool useWeighted)
{
    isOptimal = true;

    // Compute offending vars and put them into priority queue
    priority_queue<VarWithViol *, vector<VarWithViol *>, violcomp> tempQueue;
    for (int nbArcInd = 0; nbArcInd < setNBarc.size(); ++nbArcInd) {
        ArcVar *nba = setNBarc[nbArcInd];
        double rc = computeArcRedCost(nba);
        if (((nba->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((nba->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            VarWithViol *vwv = new VarWithViol();
            vwv->var = nba;
            vwv->violation = -1.0 * fabs(rc);
            tempQueue.push(vwv);
            while (tempQueue.size() > numLargestViolVarsToTrack) {
                VarWithViol *tempVWV = tempQueue.top();
                tempQueue.pop();
                delete tempVWV;
            }
        }
    }
    for (int nbEqfInd = 0; nbEqfInd < setNBeqf.size(); ++nbEqfInd) {
        EqFlowVar *nbeq = setNBeqf[nbEqfInd];
        double rc = computeEqfRedCost(nbeq);
        if (((nbeq->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((nbeq->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            VarWithViol *vwv = new VarWithViol();
            vwv->var = nbeq;
            vwv->violation = 
                useWeighted ? -1.0*fabs(rc)/(nbeq->numOrigArcs) : -1.0*fabs(rc);
            tempQueue.push(vwv);
            while (tempQueue.size() > numLargestViolVarsToTrack) {
                VarWithViol *tempVWV = tempQueue.top();
                tempQueue.pop();
                delete tempVWV;
            }
        }
    }

    vector<Variable *> finalOffenders;
    while (!tempQueue.empty()) {
        VarWithViol *tempVWV = tempQueue.top();
        tempQueue.pop();
        finalOffenders.push_back(tempVWV->var);
        delete tempVWV;
    }

    if (finalOffenders.empty()) {
        return NULL;
    } else {
        return finalOffenders[drand48() * finalOffenders.size()];
    }
}

Variable* GenNetEq::getVarBySimAnnealingArcsFirst()
{
	Variable* best;
	// Large violation means that we are (probably) earlier in the SA
	// process, so we want to prefer arcs to equal flow sets
	if (drand48() < lastViolation / baselineViolation)
		best = getVarWithLargestViolationArcsFirst();
	else
		best = getVarWithLargestViolation(true);

	// Rescale our baseline until we get close to the end of the process.
	if (lastViolation * 100 < baselineViolation && baselineViolation > 5000)
		baselineViolation = lastViolation;

	return best;
}

Variable* GenNetEq::getVarBySimAnnealingEqfsFirst()
{
	Variable* best;
	// Large violation means that we are (probably) earlier in the SA
	// process, so we want to prefer arcs to equal flow sets
	if (drand48() < lastViolation / baselineViolation)
		best = getVarWithLargestViolationEqfsFirst();
	else
		best = getVarWithLargestViolation(true);

	// Rescale our baseline until we get close to the end of the process.
//	if (lastViolation * 10 < baselineViolation && baselineViolation > 500)
//		baselineViolation = lastViolation;

	return best;
}

/*****************************************************************************/
/* Pivoting Rule Helper Functions                                            */
/*****************************************************************************/
void GenNetEq::computeArcsWithLargestViolation(double *largestViolation)
{
    for (int nbArcInd = 0; nbArcInd < setNBarc.size(); ++nbArcInd) {
        ArcVar *nba = setNBarc[nbArcInd];
        double rc = computeArcRedCost(nba);
        if (((nba->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((nba->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            checkVarForOffending(nba, fabs(rc), largestViolation);
        }
    }
    return;
}

void GenNetEq::computeEqfsWithLargestViolation(double *largestViolation, 
                                               bool useWeighted)
{
    for (int nbEqfInd = 0; nbEqfInd < setNBeqf.size(); ++nbEqfInd) {
        EqFlowVar *nbeq = setNBeqf[nbEqfInd];
        double rc = computeEqfRedCost(nbeq);
        if (((nbeq->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((nbeq->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            double frc = useWeighted ? fabs(rc) / nbeq->numOrigArcs : fabs(rc);
            checkVarForOffending(nbeq, frc, largestViolation);
        }
    }
    return;
}

inline
void GenNetEq::checkVarForOffending(Variable *var, double varViolation,
                                    double *largestViolation) 
{
    if (*largestViolation < varViolation - ROUNDOFF_TOLERANCE) {
        offendingVars.clear();
    }
    if (*largestViolation <= varViolation + ROUNDOFF_TOLERANCE) {
        *largestViolation = varViolation;
        if ((!trackLowestVarIndex) || (offendingVars.empty()) ||
            ((lastEVarInd < offendingVars.front()->varIndex) &&
             (offendingVars.front()->varIndex < var->varIndex)) ||
            ((var->varIndex < offendingVars.front()->varIndex) &&
             (offendingVars.front()->varIndex < lastEVarInd))) {
            offendingVars.push_back(var);
        } else {
            offendingVars.push_back(offendingVars.front());
            offendingVars[0] = var;
        }
    }
    return;
}

/*****************************************************************************/
/* Candidate List Pivoting Rule Functions                                    */
/*****************************************************************************/
Variable *GenNetEq::getNextCandidate(bool useWeighted)
{
    isOptimal = true;

    if ((itersSinceUpdate < 0) || (itersSinceUpdate > majorIterationLength)) {
//        printf("Rebuilding candidate list (%d)\n", itersSinceUpdate);
        updateCandidateList(useWeighted);
        itersSinceUpdate = 0;
    } else {
        ++itersSinceUpdate;
    }

    Variable *eVar = NULL;
    while (!candidateList.empty()) {
        VarWithViol *candidateVWV = candidateList.top();
        candidateList.pop();
        Variable *cVar = candidateVWV->var;
        delete candidateVWV;
        // Check if the variable still violates its optimality conditions
        bool cVarAtLB = (cVar->setLoc == L_var);
        double rc;
        if (cVar->varType == ARC_VAR) {
            rc = computeArcRedCost(dynamic_cast<ArcVar *>(cVar));
        } else { // cVar->varType == EQF_VAR
            rc = computeEqfRedCost(dynamic_cast<EqFlowVar *>(cVar));
        }
        if (((cVarAtLB) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((!cVarAtLB) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            eVar = cVar;
            break;
        }
    }

    if ((eVar == NULL) && (itersSinceUpdate > 0)) {
        // We weren't able to find anything, but we haven't updated in a 
        // while, so set flag to update and try again
        itersSinceUpdate = -1;
        return getNextCandidate();
    } else {
        // We updated this iteration but nothing violated optimality, so we 
        // return eVar, which is null, indicating that
        return eVar;
    }
}

// Running time: O(m' + np)
void GenNetEq::updateCandidateList(bool useWeighted)
{
    // First empty the candidate list and clean up memory
    clearCandidateList();

    for (int nbArcInd = 0; nbArcInd < setNBarc.size(); ++nbArcInd) {
        ArcVar *nba = setNBarc[nbArcInd];
        double rc = computeArcRedCost(nba);
        if (((nba->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((nba->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            VarWithViol *vwv = new VarWithViol();
            vwv->var = nba;
            vwv->violation = fabs(rc);
            candidateList.push(vwv);
        }
    }
    for (int nbEqfInd = 0; nbEqfInd < setNBeqf.size(); ++nbEqfInd) {
        EqFlowVar *nbeq = setNBeqf[nbEqfInd];
        double rc = computeEqfRedCost(nbeq);
        if (((nbeq->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((nbeq->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            VarWithViol *vwv = new VarWithViol();
            vwv->var = nbeq;
            vwv->violation = 
                useWeighted ? fabs(rc) / (nbeq->numOrigArcs) : fabs(rc);
            candidateList.push(vwv);
        }
    }
    return;
}

void GenNetEq::clearCandidateList()
{
    while (!candidateList.empty()) {
        VarWithViol *vwv = candidateList.top();
        candidateList.pop();
        delete vwv;
    }
    return;
}

/*****************************************************************************/
/* Candidate List 2 Pivoting Rule Functions                                  */
/*****************************************************************************/
Variable *GenNetEq::getNextCandidate2(bool useWeighted)
{
    isOptimal = true;

    if ((itersSinceUpdate < 0) || (itersSinceUpdate > majorIterationLength)) {
//        printf("Rebuilding candidate list (%d)\n", itersSinceUpdate);
        updateCandidateList2();
        itersSinceUpdate = 0;
    } else {
        ++itersSinceUpdate;
    }

    Variable *eVar = NULL;
    double largestViolation = -1.0;
    int clInd = 0;
    int eVarCLInd = -1;
    while (clInd < candidateList2.size()) {
        Variable *clVar = candidateList2[clInd];
        bool clVarAtLB = (clVar->setLoc == L_var);
        double rc;
        if (itersSinceUpdate > 0) {
            if (clVar->varType == ARC_VAR) {
                rc = computeArcRedCost(dynamic_cast<ArcVar *>(clVar));
            } else { // clVar->varType == EQF_VAR
                rc = computeEqfRedCost(dynamic_cast<EqFlowVar *>(clVar));
            }
        } else { // Reduced cost is up to date
            rc = clVar->redCost;
        }
        if (((clVarAtLB) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((!clVarAtLB) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            if ((clVar->varType == EQF_VAR) && (useWeighted)) {
                rc /= (dynamic_cast<EqFlowVar *>(clVar))->numOrigArcs;
            }
            if (largestViolation < fabs(rc)) {
                largestViolation = fabs(rc);
                eVar = clVar;
                eVarCLInd = clInd;
            }
            ++clInd; // move on to next variable
        } else { 
            // We need to remove the candidate variable from the list. This 
            // is accomplished by swapping this one with the back of the list. 
            // We do not increment the counter, though. 
            candidateList2[clInd] = candidateList2.back();
            candidateList2.pop_back();
        }
    }

    if ((eVar == NULL) && (itersSinceUpdate > 0)) {
        // We weren't able to find anything, but we haven't updated in a 
        // while, so set flag to update and try again
        itersSinceUpdate = -1;
        return getNextCandidate2();
    } else {
        // eVar may or may not be null; if it is not, we need to remove it 
        // from the candidate list, and then we can return it (regardless of 
        // whether or not it is null)
        if (eVarCLInd >= 0) {
            candidateList2[eVarCLInd] = candidateList2.back();
            candidateList2.pop_back();
        }
        return eVar;
    }
}

// Running time: O(m' + np)
void GenNetEq::updateCandidateList2()
{
    // First empty the candidate list and clean up memory
    candidateList2.clear();

    int numNBArcs = setNBarc.size();
    int numNBEqfs = setNBeqf.size();
    double arcProb = ((double) numNBArcs) / ((double) (numNBArcs + numNBEqfs));

    if (nbArcPivotInd >= numNBArcs) nbArcPivotInd = 0;
    if (nbEqfPivotInd >= numNBEqfs) nbEqfPivotInd = 0;
    int startingNBArcInd = nbArcPivotInd;
    int startingNBEqfInd = nbEqfPivotInd;
    bool exhaustedArcs = (numNBArcs == 0);
    bool exhaustedEqfs = (numNBEqfs == 0);

    // Check next batch of arcs to add to candidate list
    for (int i = 0; (i < numArcsToAlwaysCheck) && (!exhaustedArcs); ++i) {
        tryToAddArc();
        nbArcPivotInd = (nbArcPivotInd + 1) % numNBArcs;
        exhaustedArcs = (nbArcPivotInd == startingNBArcInd);
    }

    // Check next batch of equal flow sets to add to candidate list
    for (int i = 0; (i < numEqfsToAlwaysCheck) && (!exhaustedEqfs); ++i) {
        tryToAddEqf();
        nbEqfPivotInd = (nbEqfPivotInd + 1) % numNBEqfs;
        exhaustedEqfs = (nbEqfPivotInd == startingNBEqfInd);
    }

    // Now ensure that the candidate list has the appropriate size
    while ((candidateList2.size() < clMaxSize) && 
           (!exhaustedArcs) && (!exhaustedEqfs)) {
        if (drand48() < arcProb) { // Try to add an arc
            tryToAddArc();
            nbArcPivotInd = (nbArcPivotInd + 1) % numNBArcs;
            exhaustedArcs = (nbArcPivotInd == startingNBArcInd);
        } else { // Try to add an equal flow set
            tryToAddEqf();
            nbEqfPivotInd = (nbEqfPivotInd + 1) % numNBEqfs;
            exhaustedEqfs = (nbEqfPivotInd == startingNBEqfInd);
        }
    }

    if (exhaustedArcs) {
        while ((candidateList2.size() < clMaxSize) && (!exhaustedEqfs)) {
            tryToAddEqf();
            nbEqfPivotInd = (nbEqfPivotInd + 1) % numNBEqfs;
            exhaustedEqfs = (nbEqfPivotInd == startingNBEqfInd);
        }
    }

    if (exhaustedEqfs) {
        while ((candidateList2.size() < clMaxSize) && (!exhaustedArcs)) {
            tryToAddArc();
            nbArcPivotInd = (nbArcPivotInd + 1) % numNBArcs;
            exhaustedArcs = (nbArcPivotInd == startingNBArcInd);
        }
    }

    return;
}

inline 
void GenNetEq::tryToAddArc()
{
    ArcVar *nextArc = setNBarc[nbArcPivotInd];
    double rc = computeArcRedCost(nextArc);
    if (((nextArc->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
        ((nextArc->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
        isOptimal = false;
        candidateList2.push_back(nextArc);
    }
    return;
}

inline 
void GenNetEq::tryToAddEqf()
{
    EqFlowVar *nextEqf = setNBeqf[nbEqfPivotInd];
    double rc = computeEqfRedCost(nextEqf);
    if (((nextEqf->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
        ((nextEqf->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
        isOptimal = false;
        candidateList2.push_back(nextEqf);
    }
    return;
}

/*****************************************************************************/

// Running time: O(m' + np)
Variable *GenNetEq::quickSelect(bool useWeighted)
{
    isOptimal = true;
    Variable *bestVar = NULL;
    double largestViol = 0.0;

    if (nbArcPivotInd >= setNBarc.size()) nbArcPivotInd = 0;
    if (nbEqfPivotInd >= setNBeqf.size()) nbEqfPivotInd = 0;
    int nbArcPivotStart = nbArcPivotInd;
    int nbEqfPivotStart = nbEqfPivotInd;
    int arcViolations = 0;
    int eqfViolations = 0;
    bool exhaustedArcs = (setNBarc.size() == 0);
    bool exhaustedEqfs = (setNBeqf.size() == 0);

    while ((eqfViolations < 1) && !exhaustedEqfs) { 
        EqFlowVar *eqfv = setNBeqf[nbEqfPivotInd]; 
        double rc = computeEqfRedCost(eqfv);
        if (((eqfv->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((eqfv->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            ++eqfViolations;
            double frc = useWeighted ? fabs(rc) / eqfv->numOrigArcs : fabs(rc);
            if (frc > largestViol) {
                largestViol = frc;
                bestVar = eqfv;
            }
        }
        nbEqfPivotInd = (nbEqfPivotInd + 1) % setNBeqf.size();
        exhaustedEqfs = (nbEqfPivotInd == nbEqfPivotStart);
    } 

    while ((arcViolations < (numNodes)) && !exhaustedArcs) {
        ArcVar *av = setNBarc[nbArcPivotInd];
        double rc = computeArcRedCost(av);
        if (((av->setLoc == L_var) && (rc < (-1.0 * ROUNDOFF_TOLERANCE))) ||
            ((av->setLoc == U_var) && (rc > ROUNDOFF_TOLERANCE))) {
            isOptimal = false;
            ++arcViolations;
            if (fabs(rc) > largestViol) {
                largestViol = fabs(rc);
                bestVar = av;
            }
        }
        nbArcPivotInd = (nbArcPivotInd + 1) % setNBarc.size();
        exhaustedArcs = (nbArcPivotInd == nbArcPivotStart);
    }

    if ((loopIter % majorIterationLength) > 0) {
        nbArcPivotInd = nbArcPivotStart;
        nbEqfPivotInd = nbEqfPivotStart;
    }
    
    return bestVar;
}

