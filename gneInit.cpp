/*****************************************************************************/
/* Author: David Morrison                                                    */
/* Date: 2010-06-16                                                          */
/* File: gneInit.cpp                                                         */
/* Description:                                                              */
/*   Contains initialization procedures for the GenNetEq algorithm.          */
/*****************************************************************************/
#include "main.h"
#include "matrix.h"
#include "graph.h"
#include "variables.h"
#include "trees.h"
#include "genneteq.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <algorithm>

#include <vector>
using std::vector;

// Running time: O(m' + n^2 + p)
void GenNetEq::initialize(Graph *g, int initType)
{
    numNodes = (*g).getNumNodes();
    numArcs = 0;
    numEqualFlowSets = (*g).getNumEqualFlowSets();
    bigM = (*g).getBigM();
    initialType = SELFLOOPS;

    initializeNodes(g);                         // O(n)
    initializeVars(g);                          // O(m' + n^2 + p)

    //useArtificialCosts = vector<int>(vars.size(), 0); // INFEASIBLE

    formInitialBFS(g);

    phaseIiters = 0;
    phaseIIiters = 0;

    return;
}

void GenNetEq::initializeNodes(Graph *g)
{
    nodes.resize(numNodes); // This ought to ensure the right capacity 
    for (int i = 0; i < numNodes; ++i) {
        // Set up initial tree node for node i
        nodes[i].id = i;
        nodes[i].potential = 0.0;
        nodes[i].supply = (*g).supply(i);
    }
    return;
}

// Running time: O(m' + n^2 + p)
void GenNetEq::initializeVars(Graph *g)
{
    // Initialize tree nodes [O(1)]
    Tree::setStaticVariables(&nodes, &setBtI, &setBtII);

    // Initialize equal flow variables [O(p)]
    for (int r = 0; r < numEqualFlowSets; ++r) {
        EqFlowVar *efv = new EqFlowVar();
        efv->eqFlowIndex = r;
        efv->numOrigArcs = 0;
        efv->LB = 0.0;
        efv->UB = std::numeric_limits<double>::max();
        efv->actualCost = 0.0;
        efv->varType = EQF_VAR;
        for (int i = 0; i < numNodes; ++i) {
            efv->nodeVals.push_back((*g).equalFlowNodeValue(i, r));
        }
        if (r == numEqualFlowSets) efv->isArtificial = true;
        else efv->isArtificial = false;
        eqFlowVars.push_back(efv);
    }

    // Initialize arc variables [O(n^2),O(m + n)]
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNodes; ++j) {
            if ((*g).capacity(i, j) <= 0.0) { // Arc doesn't exist
                continue;
            } // Else arc exists
            if ((*g).equalFlowIndex(i, j) < 0) { // Arc not in EqFlow set
                ArcVar *av = new ArcVar();
                av->tail = i;
                av->head = j;
                av->LB = 0.0;
                av->UB = (*g).capacity(i, j);
                av->actualCost = (*g).cost(i, j);
                av->mult = (*g).multiplier(i, j);
                av->varType = ARC_VAR;
                av->isArtificial = false;
                arcVars.push_back(av);
            } else { // Arc belongs to an equal flow set
                EqFlowVar *efv = eqFlowVars[(*g).equalFlowIndex(i, j)];
                efv->actualCost += (*g).cost(i, j);
                efv->numOrigArcs += 1;
                if (efv->UB > (*g).capacity(i, j)) {
                    efv->UB = (*g).capacity(i, j);
                }
            }
        }
    }

    // Create artificial self-loops [O(n)]
    if (initialType == SELFLOOPS) {
        for (int i = 0; i < numNodes; ++i) {
            ArcVar *av = new ArcVar();
            av->tail = i;
            av->head = i;
            av->LB = 0.0;
            av->UB = std::numeric_limits<double>::max();
            av->actualCost = (*g).getBigM();
            if ((*g).supply(i) >= 0.0) {
                av->mult = 0.5;
            } else { // supply < 0
                av->mult = 2.0;
            }
            av->varType = ARC_VAR;
            av->isArtificial = true;
            arcVars.push_back(av);
        }
    }
    // We need one artifcial loop, but we don't know where it goes yet, so 
    // we'll just stick it somewhere and update it later.
    else if (initialType == INFEASIBLE) {
        ArcVar *av = new ArcVar();
        av->tail = 0;
        av->head = 0;
        av->LB = 0.0;
        av->UB = std::numeric_limits<double>::max();
        av->actualCost = (*g).getBigM();
        if ((*g).supply(0) >= 0.0) {
            av->mult = 0.5;
        } else { // supply < 0
            av->mult = 2.0;
        }
        av->varType = ARC_VAR;
        av->isArtificial = true;
        arcVars.push_back(av);
    }
    numArcs = arcVars.size(); // Not actual number of arcs, just arc variables

    // Now copy pointers to all variables into variable vector [O(m' + p)]
    int curVarIndex = 0;
    for (int aInd = 0; aInd < arcVars.size(); ++aInd) {
        ArcVar *av = arcVars[aInd];
        av->varIndex = curVarIndex;
        ++curVarIndex;
        vars.push_back(av);
    }
    for (int eqfInd = 0; eqfInd < eqFlowVars.size(); ++eqfInd) {
        EqFlowVar *eqfv = eqFlowVars[eqfInd];
        eqfv->varIndex = curVarIndex;
        ++curVarIndex;
        vars.push_back(eqfv);
    }

    return;
}

// Running time: O(m' + p)
void GenNetEq::formInitialBFS(Graph *g)
{
    // Add all of the self loops into the basis, and the remaining variables
    // into the equal flow sets
    if (initialType == SELFLOOPS)
    {
	printf("Using self-loops...\n");
        for (int avInd = 0; avInd < arcVars.size(); ++avInd) 
        {
            ArcVar *av = arcVars[avInd];
            if (av->isArtificial) 
                insertIntoB(av);
            else
                insertIntoL(av); 
        }
        for (int eqfInd = 0; eqfInd < numEqualFlowSets; ++eqfInd) {
            insertIntoL(eqFlowVars[eqfInd]);
        }
        Tree::updateTrees();
    }

    // Construct a (possibly infeasible) Type I tree, 
    else if (initialType == INFEASIBLE)
    {
        initPreds = vector<int>(numNodes, -1);
        initComps = vector<int>(numNodes, -1);
        supplyToRoot = vector<double>(numNodes, 0);
        currentCompNum = 0;
        currNumInfeas = 0;
        // Put all arcs and all equal flow sets into the lower bound set
        // initially.
        for (int i = 0; i < arcVars.size(); ++i)
            insertIntoL(arcVars[i]);
        for (int i = 0; i < numEqualFlowSets; ++i)
            insertIntoL(eqFlowVars[i]);

        int bestInMajor = (int) INFINITY;
        initCands = vector<ArcVar*>(setNBarc.begin(), setNBarc.end());
        // Continue this loop until we have a spanning tree
        while (setBarc.size() < numNodes - 1)
        {
//            printf("%d\n", setBarc.size());
            ArcVar* best = NULL;
            bool ranMajorIteration = false;

            int tries = 0;
            double bestScore = (int)INFINITY;
            while (tries < 100 || best == NULL)
            {
                int index = (int)(drand48() * initCands.size());
                double score = computeNewInfeasible(initCands[index]);
                if (score == -1)
                {
                    vector<ArcVar*>::iterator end = remove(initCands.begin(),
                            initCands.end(), initCands[index]);
                    initCands.erase(end, initCands.end());
                    continue;
                }

                else if (score < bestScore)
                {
                    bestScore = score;
                    best = initCands[index];
                }
                tries++;
            }

            // Add the best arc into the basis
            addInitialArc(best);

            //printf("(%d, %d): %0.1f\n", best->tail + 1, best->head + 1, bestInMinor);
        }

        // Now, compute the initial flow on the new-found basis structure;
        // we want to put a self-loop at the root to take care of any imbalance
        // there.
        int root;
        for (int i = 0; i < initPreds.size(); ++i)
        {
            if (initPreds[i] == -1)
            {
                root = i;
                break;
            }
        }
        // Add the self-loop into the basis, and set it to be pointing at the
        // root.
        printf("root: %d\n", root);
        for (int i = 0; i < arcVars.size(); ++i)
        {
            if (arcVars[i]->isArtificial)
            {
                arcVars[i]->tail = root;
                arcVars[i]->head = root;
                if (supplyToRoot[root] >= 0.0) 
                    arcVars[i]->mult = 0.5;
                else 
                    arcVars[i]->mult = 2.0;
                removeFromL(arcVars[i]);
                insertIntoB(arcVars[i]);
                // NOTE: Should we break here?
            }
        }

        // Reinitialize the matrix
        // NOTE: These functions are deprecated. However, this should not be 
        // needed any more. May need to check this, though...
        //clearDataStructures();
        //initializeDataStructures(g);

        // Compute the final (initial) flow from the basis structure
        Tree::updateTrees();

        computeFlows();
        int badArcs = 0;
        // NOTE: Revised this slightly; use of varIndex is safer
        for (int i = 0; i < vars.size(); ++i) {
            if (vars[i]->isArtificial && vars[i]->varType == ARC_VAR) {
                ArcVar *av = dynamic_cast<ArcVar *>(vars[i]);
                if ((av->flow < av->LB - ROUNDOFF_TOLERANCE) ||
                    (av->flow > av->UB + ROUNDOFF_TOLERANCE)) {
                    av->mult = 1.0 / av->mult;
                }
            } else if (vars[i]->flow < vars[i]->LB - ROUNDOFF_TOLERANCE) {
                badArcs++;
                useArtificialCosts[i] = -1.0;
            } else if (vars[i]->flow > vars[i]->UB + ROUNDOFF_TOLERANCE) {
                badArcs++;
                useArtificialCosts[i] = 1.0;
            }
        }
        printf("Total number of bad arcs is: %d\n", badArcs);
    }
}

// Without actually adding an arc to the tree, see how many edges it will
// make infeasible if it /was/ added to the tree.  This is tricky, because it
// needs to mimic the behaviour of addInitialArc, but not actually make any
// changes to the data structure.  There will be undefined (probably incorrect)
// behaviour if these two functions don't perform exactly the same way.  
// Unfortunately, this is almost inevitable, given the complexity of the code.
double GenNetEq::computeNewInfeasible(ArcVar* arc)
{
    int tailComp = initComps[arc->tail];
    int headComp = initComps[arc->head];

    double badness;

    // If the arc spans a single tree, this is bad
    if ((tailComp == headComp && tailComp != -1) || arc->tail == arc->head)
        return -1;

    // If neither endpoint belongs to a tree yet, it doesn't change any of the
    // rest of the tree structure, so just look to see if it's feasible or not.
    else if (tailComp == -1 && headComp == -1)
    {
        badness = score(arc, nodes[arc->tail].supply, nodes[arc->head].supply, true);
        return badness + currNumInfeas + 
            (nodes[arc->tail].supply < arc->LB || nodes[arc->tail].supply > arc->UB);
    }

    // Otherwise, at least one endpoint belongs to a tree
    else
    {
        // First, calculate all of the bad arcs before we add this arc in
        vector<int> infeasArcs(setBarc.size(), 0);
        int i = 0; 
        for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
            ArcVar *xb = setBarc[bArcInd];
            if ((xb->flow < xb->LB - ROUNDOFF_TOLERANCE) ||
                (xb->flow > xb->UB + ROUNDOFF_TOLERANCE)) {
                infeasArcs[i] = 1;
            }
            ++i;
        }

        // Walk along the path from this arc to the root of the tree
        // and compute the arcs that are made feasible or infeasible
        int curr = headComp != -1 ? arc->tail : arc->head;
        int pred = headComp != -1 ? arc->head : arc->tail;

        // Need to adjust flow down to a leaf if we're joining two trees
        double delta = 0;
        vector<double> strCopy(supplyToRoot);
        if (headComp != -1 && tailComp != -1)
        {
            vector<int> tempPred;
            double tempSupply;
            int oldCurr = arc->tail;
            int oldPred = initPreds[oldCurr];

            // Travel down to the root of the tree, adjusting the supply as
            // we go
            while (oldPred != -1)
            {
                ArcVar *avOld = findBasicAV(oldCurr, oldPred);
                if (avOld->tail == oldCurr)
                    strCopy[oldPred] -= avOld->flow * avOld->mult;
                else
                    strCopy[oldPred] += avOld->flow;
                // Keep track of the path we're walking
                tempPred.push_back(oldCurr);
                oldCurr = oldPred;
                oldPred = initPreds[oldCurr];
            }
            tempSupply = strCopy[oldCurr];

            // Walk back up, reversing the flow on all the edges, until we
            // get to the tail of the arc we're considering
            while (!tempPred.empty())
            {
                oldPred = tempPred.back();

                // Compute the flow, assuming that no flow is going over the
                // edge (this is the false value).  This is equivalent to 
                // setting the flow to 0 in the actual flow calculation routine
                computeInitFlowInfeasible(oldCurr, oldPred, tempSupply, 
                        false, strCopy, infeasArcs);
                tempSupply = strCopy[oldPred];
                oldCurr = oldPred;
                tempPred.pop_back();
            }
            delta = tempSupply;

            badness = score(arc, strCopy[arc->tail], strCopy[arc->head], true);

            // Arc is not in the basis, so we need to check this one by hand
            if (delta < arc->LB || delta > arc->UB)
                badness += 1.0;
            
            delta *= arc->mult;
            delta += supplyToRoot[pred];
            curr = arc->head;
            pred = initPreds[curr];
        }

        else if (tailComp == -1)
        {
            // This is what the flow on the arc would be (need to check by hand,
            // since the arc is not actually in the basis yet)
            badness = score(arc, strCopy[arc->tail], strCopy[arc->head], true);
            double oldSupply = strCopy[pred];
            delta = nodes[arc->tail].supply;
            if (delta < arc->LB || delta > arc->UB)
                badness += 1.0;

            // This is what the supply at the next node up the tree would be
            strCopy[pred] += delta * arc->mult;
            delta = strCopy[pred] - oldSupply;
            
            // Move along the path to the root
            curr = arc->head;
            pred = initPreds[curr];
        }
        else if (headComp == -1)
        {
            // This is what the flow on the arc would be (need to check by hand,
            // since the arc is not actually in the basis yet)
            badness = score(arc, strCopy[arc->head], strCopy[arc->tail], false);
            double oldSupply = strCopy[pred];
            delta = nodes[arc->head].supply;
            delta /= -(arc->mult);
            if (delta < arc->LB || delta > arc->UB)
                badness += 1.0;

            // This is what the supply at the next node up the tree would be
            strCopy[pred] += delta / arc->mult;
            delta = strCopy[pred] - oldSupply;

            // Move along the path to the root
            curr = arc->tail;
            pred = initPreds[curr];
        }

        // Compute new flow along the path to the root
        while (pred != -1)
        {
            double oldSupply = strCopy[pred];
            computeInitFlowInfeasible(curr, pred, delta, true,
                    strCopy, infeasArcs);
            delta = strCopy[pred] - oldSupply;
            curr = pred;
            pred = initPreds[curr];
        }

        // Count up the total number of infeasible arcs
        for (int i = 0; i < infeasArcs.size(); ++i)
            badness += infeasArcs[i];

        return badness;
    }
}

// Add an arc to the initial basis structure, and adjust the supply/flow at all
// affected nodes and arcs, as well as maintain the necessary tree structure.
void GenNetEq::addInitialArc(ArcVar* arc)
{
    removeFromL(arc);
    insertIntoB(arc);
    int tailComp = initComps[arc->tail];
    int headComp = initComps[arc->head];
    int curr, pred;
    double delta;

    // If the arc spans a single tree, this is bad
    if ((tailComp == headComp && tailComp != -1) || arc->tail == arc->head)
    {
        fprintf(stderr, "ERROR: Trying to add bad arc.  Something is wrong.\n");
        return;
    }

    // If neither endpoint belongs to a tree yet, just add it
    if (tailComp == -1 && headComp == -1)
    {
        // Create a new tree
        initComps[arc->tail] = currentCompNum;
        initComps[arc->head] = currentCompNum;
        initPreds[arc->tail] = arc->head;
        currentCompNum++;
        curr = arc->tail;
        pred = arc->head;
        delta = nodes[curr].supply;
        supplyToRoot[curr] = delta;
        supplyToRoot[pred] = nodes[pred].supply;
    }

    // If both endpoints belong to a tree, we need to merge the trees
    else if (tailComp != -1 && headComp != -1)
    {
        // Merge the trees together
        int old = initComps[arc->tail];
        initComps[arc->tail] = initComps[arc->head];
        for (int i = 0; i < initComps.size(); ++i)
            if (initComps[i] == old)
                initComps[i] = initComps[arc->head];

        int oldPred = initPreds[arc->tail];
        int oldCurr = arc->tail;

        // Reverse all the predecessor labels
        while (oldPred != -1)
        {
            // While we're doing this, remove any effect that the current
            // arc has on the predecessor supply -- however, we want to ensure
            // that the supply effects of all other nodes in the tree remain
            // constant.  Note that we don't do anything with the current arcs
            // yet; everything will get touched except arc->tail, and will get
            // modified appropriately on the way back up.
            ArcVar *avOld = findBasicAV(oldCurr, oldPred);
            if (avOld->tail == oldCurr)
                supplyToRoot[oldPred] -= avOld->flow * avOld->mult;
            else
                supplyToRoot[oldPred] += avOld->flow;

            // Because computeInitFlowOnArc is additive with respect to the
            // flow, we need to clear the flows along this path
            avOld->flow = 0.0;

            // Swap the predecessor labels
            int temp = initPreds[oldPred];
            initPreds[oldPred] = oldCurr;
            oldCurr = oldPred;
            oldPred = temp;
        }

        // Recompute the flow along this path -- delta is left over from the
        // previous while loop, and should be set to the supply value of the
        // old root.
        while (oldCurr != arc->tail)
        {
            oldPred = initPreds[oldCurr];
            delta = supplyToRoot[oldCurr];
            computeInitFlowOnArc(oldCurr, oldPred, delta);
            oldCurr = oldPred;
        }
        initPreds[arc->tail] = arc->head;
        curr = arc->tail;
        pred = arc->head;

        // NOT (*supply)[curr]
        delta = supplyToRoot[curr];
    }

    // If only one endpoint belongs to a tree, just add the other to be the same
    else if (tailComp == -1)
    {
        initComps[arc->tail] = initComps[arc->head];
        initPreds[arc->tail] = arc->head;
        curr = arc->tail;
        pred = arc->head;
        delta = nodes[curr].supply;
        supplyToRoot[curr] = delta;
    }
    else if (headComp == -1)
    {
        initComps[arc->head] = initComps[arc->tail];
        initPreds[arc->head] = arc->tail;
        curr = arc->head;
        pred = arc->tail;
        delta = nodes[curr].supply;
        supplyToRoot[curr] = delta;
    }

    // The only flow in the graph that has changed is the flow along the path
    // from the current arc to the root of the tree.  So run through these arcs
    // and update.
    while (pred != -1)
    {
        double oldSupply = supplyToRoot[pred];
        computeInitFlowOnArc(curr, pred, delta);
        delta = supplyToRoot[pred] - oldSupply;
        curr = pred;
        pred = initPreds[curr];
    }

    currNumInfeas = computeNumInfeasible();
}

double GenNetEq::score(ArcVar* arc, double su, double sv, bool forward)
{
    double mu = arc->mult;
    double cap = arc->UB;

    double supplyU = forward ? su : -su/mu;
    double supplyV = forward ? -sv/mu : sv;

    // We only need to check at a few critical points: the lower bound, the
    // upper bound, and when the flow at the head/tail nodes is completely
    // satisfied.
    vector<double> critPoints;
    critPoints.push_back(0);
    critPoints.push_back(cap);
    if (supplyU >= 0 && supplyU <= cap)
        critPoints.push_back(supplyU);

    if (supplyV >= 0 && supplyV <= cap)
        critPoints.push_back(supplyV);

    double maxSatisfied = -INFINITY;
    double flowAtMax = 0.0;
    for (int i = 0; i < critPoints.size(); ++i)
    {
        double atU, atV;
        if (forward)
        {
            atU = supplyU > 0 ? supplyU - fabs(supplyU - critPoints[i]) :
                supplyU - fabs(critPoints[i]);
            atV = supplyV * mu > 0 ? supplyV * mu - fabs(mu * critPoints[i]) :
                supplyV * mu - fabs(supplyV * mu - mu * critPoints[i]);
        }
        else
        {
            atU = supplyU * mu > 0 ? supplyU * mu - fabs(mu * critPoints[i]) :
                supplyU * mu - fabs(supplyU * mu - mu * critPoints[i]);
            atV = supplyV > 0 ? supplyV - fabs(supplyV - critPoints[i]) :
                supplyV - fabs(critPoints[i]);
        }

        if (atU + atV > maxSatisfied)
        {
            maxSatisfied = atU + atV;
            flowAtMax = critPoints[i];
        }
    }

    return exp(-maxSatisfied / 5.0);
}

int GenNetEq::computeNumInfeasible()
{
    // Find out how many arcs in the tree are infeasible with 
    // this choice of new arc
    int badArcs = 0;
    for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
        ArcVar *xb = setBarc[bArcInd];
        if ((xb->flow < xb->LB - ROUNDOFF_TOLERANCE) ||
            (xb->flow > xb->UB + ROUNDOFF_TOLERANCE)) {
            badArcs++;
        }
    }
    return badArcs;
}

// Adjust the flow on an arc, assuming that the supply at curr changed by
// delta.  Update the supply at pred with the new value.
void GenNetEq::computeInitFlowOnArc(int curr, int pred, double delta)
{
    ArcVar* a = findBasicAV(curr, pred);
    int tail = a->tail;
    int head = a->head;
    double mult = a->mult;

    a->flow += (curr == tail) ? delta : -delta / mult;
    supplyToRoot[pred] += (curr == tail) ? mult * delta : delta / mult;
}

// See if the flow on arc would become infeasible if the supply at curr changed
// by delta.  Update the supply at the predecessor to what it would be in this
// case; if we're having to reverse the direction of a path, we want to reset
// the flow to 0 before computing, so set useFlow to false.
void GenNetEq::computeInitFlowInfeasible(int curr, int pred, double delta, 
        bool useFlow, vector<double>& strCopy, vector<int>& infeasArcs)
{
    int basisInd = 0;
    ArcVar *arc = NULL;
    for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
        ArcVar *xb = setBarc[bArcInd];
        if ((xb->tail == curr && xb->head == pred) || 
            (xb->tail == pred && xb->head == curr)) {
            arc = xb;
            break;
        }
        ++basisInd;
    }

    double flow = useFlow ? arc->flow : 0.0;
    
    if (curr == arc->tail)
    {
        if (flow + delta < arc->LB || flow + delta > arc->UB)
            infeasArcs[basisInd] = 1;
        else
            infeasArcs[basisInd] = 0;
        strCopy[pred] += arc->mult * delta;
    }
    else
    {
        delta /= -(arc->mult);
        if (flow + delta < arc->LB || flow + delta > arc->UB)
            infeasArcs[basisInd] = 1;
        else
            infeasArcs[basisInd] = 0;
        strCopy[pred] += delta / arc->mult;
    }
}

// Return the pointer to the arc variable starting at u and ending at v
ArcVar* GenNetEq::findAV(int u, int v)
{
    for (int i = 0; i < arcVars.size(); ++i)
        if (arcVars[i]->tail == u && arcVars[i]->head == v)
            return arcVars[i];

    return NULL;
}

// Return the pointer to the basic arc variable with one endpoint at u
// and the other at v (not necessarily in that order)
ArcVar* GenNetEq::findBasicAV(int u, int v)
{
    for (int bArcInd = 0; bArcInd < setBarc.size(); ++bArcInd) {
        ArcVar *xb = setBarc[bArcInd];
        if ((xb->tail == u && xb->head == v) || 
            (xb->tail == v && xb->head == u)) {
            return xb;
        }
    }
    return NULL;
}

