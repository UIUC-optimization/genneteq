/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-06-16                                                          */
/* File: gnePivot.cpp                                                        */
/* Description:                                                              */
/*   Contains implementation details for the pivoting functions used by the  */
/*   GenNetEq algorithm.                                                     */
/*****************************************************************************/
#include "main.h"
#include "variables.h"
#include "trees.h"
#include "genneteq.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>

#include <vector>
using std::vector;
#include <list>
using std::list;

/*****************************************************************************/
/* Pivoting Protected Functions                                              */
/*****************************************************************************/
// Running time: O(ns + g(s))
void GenNetEq::pivot() 
{
    // Compute yFlow [O(ns + g(s))]
    if ((eVar->varType == ARC_VAR) && 
        (nodes[dynamic_cast<ArcVar *>(eVar)->tail].memberOfTreeType==TYPE_I)&&
        (nodes[dynamic_cast<ArcVar *>(eVar)->head].memberOfTreeType==TYPE_I)){
        // If the entering variable is an arc that crosses two type I trees, 
        // then the equal flow sets will not change and we can save computing 
        // time by ignoring them
        pivotTI(dynamic_cast<ArcVar *>(eVar));
    } else {
        // Otherwise the entering variable is either an arc that involves at 
        // least one type II tree or an equal flow set, in which case the 
        // basic equal flow sets will be needed in computing the new flow
        pivotTII();
    }
    if (lVar == NULL) {
        return;
    }

    // DEBUG:
//    printf("Iter %6d; EVar: %d; LVar: %d\n", loopIter, eVar->varIndex, 
//                                                       lVar->varIndex);

//    if ((initialType == INFEASIBLE) && presolve) updateArtificialCosts(eVar); INFEASIBLE
    updateBasis();

    // Check for cycling, and update pivoting stats [O(|C|)]
    checkForCycling();
    updatePivotingStats();

    return;
}

// Running time: O(n)
void GenNetEq::pivotTI(ArcVar *eav)
{
    TreeTypeI *tU = setBtI[nodes[eav->tail].memberOfTreeID];
    TreeTypeI *tV = setBtI[nodes[eav->head].memberOfTreeID]; 

    // Compute blocking variables [O(n)]
    delta = std::numeric_limits<double>::max(); // Will clear blocking vars
    computeFlowsPivotingTI(eav, tU, tV);

    // Identify leaving variable
    identifyLeavingVar(); 
    if (lVar == NULL) {
        return;
    }

    // Update flow based on delta [O(n)] 
    vector<ArcVar *> *tArcs = tU->getArcsInTree();
    for (int tArcInd = 0; tArcInd < (*tArcs).size(); ++tArcInd) {
        ArcVar *av = (*tArcs)[tArcInd];
        av->flow += (delta * (av->deltaFlow));
        flowCost += av->cost * (delta * (av->deltaFlow));
        av->deltaFlow = 0.0;
//        checkVarFlow(av); // Error-checking
    }
    ArcVar *exav = tU->getExtraArc();
    exav->flow += (delta * (exav->deltaFlow));
    flowCost += exav->cost * (delta * (exav->deltaFlow));
    exav->deltaFlow = 0.0;
//    checkVarFlow(exav); // Error-checking

    if (tU->getTreeID() != tV->getTreeID()) {
        tArcs = tV->getArcsInTree();
        for (int tArcInd = 0; tArcInd < (*tArcs).size(); ++tArcInd) {
            ArcVar *av = (*tArcs)[tArcInd];
            av->flow += (delta * (av->deltaFlow));
            flowCost += av->cost * (delta * (av->deltaFlow));
            av->deltaFlow = 0.0;
//            checkVarFlow(av); // Error-checking
        }
        exav = tV->getExtraArc();
        exav->flow += (delta * (exav->deltaFlow));
        flowCost += exav->cost * (delta * (exav->deltaFlow));
        exav->deltaFlow = 0.0;
//        checkVarFlow(exav); // Error-checking
    }
    // Handle entering variable's flow
    if (eav->setLoc == L_var) { // eVar at LB
        eav->flow = delta;
        flowCost += eav->cost * delta;
    } else { // eVar at UB
        eav->flow = eav->UB - delta;
        flowCost -= eav->cost * delta;
    }
//    checkVarFlow(eav); // Error-checking

    return;
}

// Running time: O(ns + g(s))
void GenNetEq::pivotTII()
{
    // Compute blocking variables [O(n)]
    delta = std::numeric_limits<double>::max(); // Will clear blocking vars
    computeFlowsPivotingTII(eVar); // [O(ns + g(s))]

    // Identify leaving variable
    identifyLeavingVar();
    if (lVar == NULL) {
        return;
    }

    // Update flow based on delta [O(n)] 
    for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
        ArcVar *av = setBarc[bArcInd];
        av->flow += (delta * (av->deltaFlow));
        flowCost += av->cost * (delta * (av->deltaFlow));
        av->deltaFlow = 0.0;
//        checkVarFlow(av); // Error-checking
    }
    for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
        EqFlowVar *eqfv = setBeqf[bEqfInd];
        eqfv->flow += (delta * (eqfv->deltaFlow));
        flowCost += eqfv->cost * (delta * (eqfv->deltaFlow));
        eqfv->deltaFlow = 0.0;
//        checkVarFlow(eqfv); // Error-checking
    }
    // Handle entering variable's flow
    if (eVar->setLoc == L_var) { // eVar at LB
        eVar->flow = delta;
        flowCost += eVar->cost * delta;
    } else { // eVar at UB
        eVar->flow = eVar->UB - delta;
        flowCost -= eVar->cost * delta;
    }
//    checkVarFlow(eVar); // Error-checking

    return;
}

/*****************************************************************************/
/* Helper functions                                                          */
/*****************************************************************************/
// Running time: O(1)
void GenNetEq::identifyLeavingVar()
{
    if (eVar->UB <= delta + ROUNDOFF_TOLERANCE) { // Entering var is blocking
        delta = eVar->UB;
        lVar = eVar;
    } else { // Leaving variable is one of the basic blocking variables
        if (blockingVars.empty()) {
            isUnbounded = true;
            lVar = NULL;
        }
        if (selectRandomBlockingVar) {
            lVar = blockingVars[drand48() * blockingVars.size()];
        } else {
            lVar = blockingVars.front(); // Take the first blocking variable
        }
    }

    // Update stats based on pivot
    if (delta <= ROUNDOFF_TOLERANCE) {
        ++totalDegenPivots;
        ++consecutiveDegenPivots;
    } else {
        consecutiveDegenPivots = 0;
    }

    // DEBUG:
//    if (lVar != NULL) {
//        double newLVarFlow = lVar->flow + delta * lVar->deltaFlow;
//        if ((fabs(newLVarFlow - lVar->LB) > ROUNDOFF_TOLERANCE) && 
//            (fabs(newLVarFlow - lVar->UB) > ROUNDOFF_TOLERANCE)) {
//            // Var is not close to either LB or UB, so we have roundoff error
//            double deltaLVar = computeDeltaVar(lVar);
//            printf("ERROR (%d): Computed delta does not take lVar %d to a bound\n",
//                loopIter, lVar->varIndex);
//            printf("%.16lf + %.16lf * %.16lf = %.16lf\n", lVar->flow, delta, 
//                lVar->deltaFlow, newLVarFlow);
//            printf("lVar Delta := %.16lf, Delta := %.16lf, Diff := %.16lf\n", 
//                deltaLVar, delta, (deltaLVar - delta));
//            exit(-1);
//        }
//    }

    return;
}

/*****************************************************************************/
/* Standard updating functions                                               */
/*****************************************************************************/
// Running time: O(n)
void GenNetEq::updateArtificialCosts() 
{
    // Handle basic arcs [O(n)]
    for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
        ArcVar *av = setBarc[bArcInd];
        if (((useArtificialCosts[av->varIndex] < 0) && 
             (av->flow > av->LB - ROUNDOFF_TOLERANCE)) ||
            ((useArtificialCosts[av->varIndex] > 0) && 
             (av->flow < av->UB + ROUNDOFF_TOLERANCE))) {
            useArtificialCosts[av->varIndex] = 0.0;
            av->cost = 0.0;
        } 
        checkVarFlow(av); // Error-checking
    }
    // Handle basic equal flow sets [O(n)]
    for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
        EqFlowVar *eqfv = setBeqf[bEqfInd];
        if (((useArtificialCosts[eqfv->varIndex] < 0) && 
             (eqfv->flow > eqfv->LB - ROUNDOFF_TOLERANCE)) ||
            ((useArtificialCosts[eqfv->varIndex] > 0) && 
             (eqfv->flow < eqfv->UB + ROUNDOFF_TOLERANCE))) {
            useArtificialCosts[eqfv->varIndex] = 0.0;
            eqfv->cost = 0.0;
        }
        checkVarFlow(eqfv); // Error-checking 
    }
    // Handle entering variable
    if (((useArtificialCosts[eVar->varIndex] < 0) && 
         (eVar->flow > eVar->LB - ROUNDOFF_TOLERANCE)) ||
        ((useArtificialCosts[eVar->varIndex] > 0) && 
         (eVar->flow < eVar->UB + ROUNDOFF_TOLERANCE))) {
        useArtificialCosts[eVar->varIndex] = 0.0;
        eVar->cost = 0.0;
    } 
    checkVarFlow(eVar); // Error-checking 
    return;
}

// Running time: O(n)
void GenNetEq::updateBasis()
{
    if (eVar->varIndex == lVar->varIndex) { // O(1)
        // When the entering variable is the same as the leaving variable, 
        // updating the basis is relatively straightforward
        if (eVar->setLoc == L_var) { // eVar at LB 
            //eVar->setLoc = U_var;
            removeFromL(eVar);
            insertIntoU(eVar);
        } else { // eVar at UB
            //eVar->setLoc = L_var;
            removeFromU(eVar);
            insertIntoL(eVar);
        }
    } else { // O(n)
        // Otherwise the entering and leaving variables differ, so we need to 
        // be more careful in how we update things
        if (eVar->setLoc == L_var) { // eVar at LB
            removeFromL(eVar);
        } else { // eVar at UB
            removeFromU(eVar);
        }

        // Swap in basis
        removeFromB(lVar);
        insertIntoB(eVar);

        // Need to determine where to put the leaving variable. One way to do 
        // this is to check if it's flow is less than half of its upper bound.
        if (lVar->flow > (0.5 * (lVar->UB))) { // lVar to UB
            insertIntoU(lVar);
        } else { // lVar to LB
            insertIntoL(lVar);
        }

        Tree::updateTrees(); // [O(n)]
        // Error-checking to ensure that trees are being updated properly
        //Tree::checkTrees(); // [O(# trees)]
    }
    // DEBUG: Now, ensure that flow on lVar is at UB or LB
    if (lVar->setLoc == L_var) {
        if (fabs(lVar->flow - lVar->LB) > ROUNDOFF_TOLERANCE) {
            printf("ERROR (%d): Correcting roundoff in flow for L variable\n", 
                loopIter);
            lVar->flow = lVar->LB;
        }
    } else { // lVar->setLoc == U_var
        if (fabs(lVar->flow - lVar->UB) > ROUNDOFF_TOLERANCE) {
            printf("ERROR (%d): Correcting roundoff in flow for U variable\n", 
                loopIter);
            lVar->flow = lVar->UB;
        }
    }

    return;
}

// Running time: O(|C|) (|C| is the size of the cycle list)
void GenNetEq::checkForCycling()
{
    // Keep the 10 most recent states to check against for cycling.  Most
    // observed cycling only takes place within 2-4 states, so this should
    // be sufficient.  Bump this up later if need be.
    if (cycleList.size() >= 10) {
        PivotVars *old = cycleList.back();
        cycleList.pop_back();
        delete old;
    }
    // Store the most recent pivot for cycle-checking
    PivotVars *curr = new PivotVars();
    curr->eVarInd = eVar->varIndex;
    curr->lVarInd = lVar->varIndex;
    cycleList.push_front(curr);

    cyclingDetected = false;
    for (list<PivotVars *>::iterator i = ++cycleList.begin(); 
            i != cycleList.end(); ++i) {
        if ((*cycleList.begin())->eVarInd == (*i)->eVarInd &&
            (*cycleList.begin())->lVarInd == (*i)->lVarInd) {
            cyclingDetected = true;
            break;
        }
    }
    // Increment stats and check for consecutive cycling. 
    if (cyclingDetected) {
        ++totalCyclesDetected;
        ++consecutiveCycles;
//        printf("Cycle detected!\n");
        if (consecutiveCycles > 5) {
            printf("Stuck in cycling loop.  Retry with different pivots.\n");
            exit(-1);
        } 
    } else { // We've broken a string of consecutive cycles
//        if (numConsecutiveCycles > 1) {
//            printf("Consecutive cycling stopped\n");
//        }
        consecutiveCycles = 0;
    }

    lastEVarInd = eVar->varIndex;

    return;
}

// Running time: O(1)
void GenNetEq::updatePivotingStats()
{
    // Stats-tracking
    totalPivots++;

    if (eVar->varType == ARC_VAR) {
        int u = dynamic_cast<ArcVar *>(eVar)->tail;
        int v = dynamic_cast<ArcVar *>(eVar)->head;
        if (lVar->varType == ARC_VAR) {
            arcArcPivots++;
        } else if (lVar->varType == EQF_VAR) {
            arcEFSPivots++;
        }

        if (nodes[u].memberOfTreeType == nodes[v].memberOfTreeType) {
            if (nodes[u].memberOfTreeID == nodes[v].memberOfTreeID &&
                        nodes[u].memberOfTreeType == TYPE_I) {
                avgTIsize += setBtI[nodes[u].memberOfTreeID]->getSize();
                sameTIPivots++;
            } else if (nodes[u].memberOfTreeID == nodes[v].memberOfTreeID &&
                        nodes[u].memberOfTreeType == TYPE_II) {
                avgTIsize += setBtII[nodes[u].memberOfTreeID]->getSize(); 
                sameTIIPivots++;
            } else if (nodes[u].memberOfTreeType == TYPE_I) {
                tItIPivots++;
            } else if (nodes[u].memberOfTreeType == TYPE_II) {
                tIItIIPivots++;
            }
        } else { 
            tItIIPivots++;
        }
    }
    if (eVar->varType == EQF_VAR) {
        EFSPivots++;
        if (lVar->varType == ARC_VAR) {
            EFSArcPivots++;
        } else if (lVar->varType == EQF_VAR) {
            EFSEFSPivots++;
        }
    }
    return;
}

/*****************************************************************************/
/* Used for steepest descent pivoting                                        */
/*****************************************************************************/
// Running time: O(ns + g(s))
// NOTE: Currently out of date 
double GenNetEq::computeMinDelta(Variable *var) 
{
    // Used to compute allowable movement when pivoting on given variable
    // Compute yFlow [O(ns + g(s))]
    if ((var->varType == ARC_VAR) && 
        (nodes[dynamic_cast<ArcVar *>(var)->tail].memberOfTreeType==TYPE_I)&&
        (nodes[dynamic_cast<ArcVar *>(var)->head].memberOfTreeType==TYPE_I)){
        // If the entering variable is an arc that crosses two type I trees, 
        // then the equal flow sets will not change and we can save computing 
        // time by ignoring them
        int uID = nodes[dynamic_cast<ArcVar *>(var)->tail].memberOfTreeID;
        int vID = nodes[dynamic_cast<ArcVar *>(var)->head].memberOfTreeID;
        computeFlowsPivotingTI(dynamic_cast<ArcVar *>(var), setBtI[uID], 
                                                             setBtI[vID]);
    } else {
        // Otherwise the entering variable is either an arc that involves at 
        // least one type II tree or an equal flow set, in which case the 
        // basic equal flow sets will be needed in computing the new flow
        computeFlowsPivotingTII(var);
    }

    // Compute minimum delta for this entering variable. By default the delta 
    // value is constrained to be no larger than the variable's upper bound. 
    double delta = var->UB;

    // Handle basic arcs [O(n)]
    for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
        ArcVar *av = setBarc[bArcInd];
        double deltaArc = computeDeltaVar(av);
        if (delta > deltaArc) {
            delta = deltaArc;
        }
        av->deltaFlow = 0.0;
    }
    // Handle basic equal flow sets [O(n)]
    for (int bEqfInd = 0; bEqfInd < setBeqf.size(); ++bEqfInd) {
        EqFlowVar *eqfv = setBeqf[bEqfInd];
        double deltaEqf = computeDeltaVar(eqfv);
        if (delta > deltaEqf) {
            delta = deltaEqf;
        }
        eqfv->deltaFlow = 0.0;
    }

    return delta;
}

/*****************************************************************************/
/* Miscellaneous Pivoting Functions                                          */
/*****************************************************************************/
// Running time: O(n^2)
void GenNetEq::convertBasicArcsToTrees()
{
    // Clear out any previous information in trees [O(n)]
    for (int i = 0; i < numNodes; ++i) {
        nodes[i].memberOfTreeID = NULL_TREE_ID;
        nodes[i].memberOfTreeType = NULL_TREE_TYPE;

        nodes[i].depth = -1;
        nodes[i].thread = -1;
        nodes[i].pred = -1;
        nodes[i].predArc = NULL;
    }

    Tree::setStaticVariables(&nodes, &setBtI, &setBtII);
    // Now add basic arcs into trees, building up forest arc by arc [O(n^2)]
    for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
        Tree::insertIntoTrees(setBarc[bArcInd]);
    }
    // Now need to check to see if we have any stranded nodes // [O(n)]
    Tree::checkForStrandedNodes();
    // Finally, update tree threads and traversal order [O(n)]
    Tree::updateTrees();
    return;
}

/*****************************************************************************/
/* Pivoting Functions for Updating Basis Structure as Vectors (pop_back)     */
/*****************************************************************************/
// Running time: O(n)
void GenNetEq::insertIntoB(Variable *var)
{
    // Variable should be added to B
    if (var->varType == ARC_VAR) {
        setBarc.push_back(dynamic_cast<ArcVar *>(var));
        var->setLoc = B_var;
        var->setInd = setBarc.size() - 1;
        Tree::insertIntoTrees(dynamic_cast<ArcVar *>(var));
    } else { // var->varType == EQF_VAR
        setBeqf.push_back(dynamic_cast<EqFlowVar *>(var));
        var->setLoc = B_var;
        var->setInd = setBeqf.size() - 1;
    }
    return;
}

// Running time: O(n)
void GenNetEq::removeFromB(Variable *var)
{
    // Variable should be removed from B
    if (var->varType == ARC_VAR) {
        int prevInd = var->setInd;
        ArcVar *varToMove = setBarc.back();
        setBarc[prevInd] = varToMove;
        varToMove->setInd = prevInd;
        setBarc.pop_back();
        Tree::removeFromTrees(dynamic_cast<ArcVar *>(var));
    } else { // var->varType == EQF_VAR
        int prevInd = var->setInd;
        EqFlowVar *varToMove = setBeqf.back();
        setBeqf[prevInd] = varToMove;
        varToMove->setInd = prevInd;
        setBeqf.pop_back();
    }
    var->setLoc = NIL;
    var->setInd = -1;
    return; 
}

// Running time: O(1)
void GenNetEq::insertIntoL(Variable *var)
{ 
    // Variable should be added to L
    if (var->varType == ARC_VAR) {
        setNBarc.push_back(dynamic_cast<ArcVar *>(var));
        var->setLoc = L_var;
        var->setInd = setNBarc.size() - 1;
    } else { // var->varType == EQF_VAR
        setNBeqf.push_back(dynamic_cast<EqFlowVar *>(var));
        var->setLoc = L_var;
        var->setInd = setNBeqf.size() - 1;
    }
    return;
}

// Running time: O(1)
void GenNetEq::removeFromL(Variable *var)
{
    // Variable should be removed from L
    if (var->varType == ARC_VAR) {
        int prevInd = var->setInd;
        ArcVar *varToMove = setNBarc.back();
        setNBarc[prevInd] = varToMove;
        varToMove->setInd = prevInd;
        setNBarc.pop_back();
    } else { // var->varType == EQF_VAR
        int prevInd = var->setInd;
        EqFlowVar *varToMove = setNBeqf.back();
        setNBeqf[prevInd] = varToMove;
        varToMove->setInd = prevInd;
        setNBeqf.pop_back();
    }
    var->setLoc = NIL;
    var->setInd = -1;
    return;
}

// Running time: O(1)
void GenNetEq::insertIntoU(Variable *var)
{
    // Variable should be added to U
    if (var->varType == ARC_VAR) {
        setNBarc.push_back(dynamic_cast<ArcVar *>(var));
        var->setLoc = U_var;
        var->setInd = setNBarc.size() - 1;
    } else { // var->varType == EQF_VAR
        setNBeqf.push_back(dynamic_cast<EqFlowVar *>(var));
        var->setLoc = U_var;
        var->setInd = setNBeqf.size() - 1;
    }
    return;
}

// Running time: O(1)
void GenNetEq::removeFromU(Variable *var)
{
    // Variable should be removed from U
    if (var->varType == ARC_VAR) {
        int prevInd = var->setInd;
        ArcVar *varToMove = setNBarc.back();
        setNBarc[prevInd] = varToMove;
        varToMove->setInd = prevInd;
        setNBarc.pop_back();
    } else { // var->varType == EQF_VAR
        int prevInd = var->setInd;
        EqFlowVar *varToMove = setNBeqf.back();
        setNBeqf[prevInd] = varToMove;
        varToMove->setInd = prevInd;
        setNBeqf.pop_back();
    }
    var->setLoc = NIL;
    var->setInd = -1;
    return;
}

