/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-06-16                                                          */
/* File: graph.cpp                                                           */
/* Description:                                                              */
/*   Contains implementation details for a custom graph class.               */
/*****************************************************************************/
#include "main.h"
#include "matrix.h"
#include "graph.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <limits>
#include <fstream>

Graph::Graph()
{
    // Default constructor, do nothing
}

Graph::~Graph() 
{
    for (int i = 0; i < numNodes; ++i) {
        delete nodes[i];
    }
}

void Graph::initialize(Config *conf) 
{
    numNodes = 0;
    numArcs = 0;
    numEqualFlowSets = 0;
    totalSupply = 0.0;
    maxArcCost = 0.0;
    maxArcCapacity = 0.0;

    readColFile(conf->inFile);
//    readGSFile(conf->inFile);

    bigM = numArcs * maxArcCost * maxArcCapacity;

	/*if (conf->initType != SELFLOOPS)
	{
		computeBestSpanForest(conf->initType);
		if (conf->initType == EQSELFLOOPS)
			addSelfLoops();
	}*/

    return;
}

void Graph::addSelfLoop(int i)
{
	printf("Adding loop (%d %d)\n", i + 1);
	// Add a self loop here
	if (nodes[i]->supply > 1)
		arcMultipliers.set(i, i, 0.5);
	else
		arcMultipliers.set(i, i, 2);
	eqFlowSetIndices.set(i, i, -1);
	arcCosts.set(i, i, bigM);
	arcCapacities.set(i, i, INFINITY);
}

void Graph::print() const 
{
    for (int i = 0; i < numNodes; ++i) {
        printf("Node %d: Index:=%d Supply:=%lf [", 
                nodes[i]->id, nodes[i]->index, nodes[i]->supply);
        for (int r = 0; r < numEqualFlowSets; ++r) {
            printf("d_%d(%d):=%lf,", r + EQF_IND_OFFSET, 
                    i + NODE_IND_OFFSET, eqFlowNodeValues.get(i, r));
        }
        printf("]\n");
    }
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNodes; ++j) {
            if ((i == j) || (arcCapacities.get(i, j) <= 0)) {
                continue;
            }
            printf("(%d -> %d): Cost:=%lf Cap:=%lf Mult:=%lf " 
                   "EqFlowSet:=%d\n", i + NODE_IND_OFFSET, j + NODE_IND_OFFSET,
                    arcCosts.get(i, j), arcCapacities.get(i, j), 
                    arcMultipliers.get(i, j), 
                    eqFlowSetIndices.get(i, j) + EQF_IND_OFFSET);
        }
    }
    return;
}

void Graph::exportColFile(const char *colFileOut) const
{
    FILE *cfo = fopen(colFileOut, "w");
    if (cfo == NULL) {
        printf("Unable to open file %s\n", colFileOut);
        exit(-1);
    }
    fprintf(cfo, "c   n1  n2      LB  UB  COST    MULT    EQFLW\n");
    fprintf(cfo, "p efpgn %d %d %d\n", numNodes, numArcs, numEqualFlowSets);
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNodes; ++j) {
            if (capacity(i, j) <= 0.0) { // Arc doesn't exist
                continue; 
            } // else arc exists
            fprintf(cfo, "a %d %d %lf %lf %lf %lf %d\n", 
                i + NODE_IND_OFFSET, j + NODE_IND_OFFSET, 0.0, 
                capacity(i, j), cost(i, j), multiplier(i, j), 
                equalFlowIndex(i, j) + EQF_IND_OFFSET);
        }
    }
    for (int i = 0; i < numNodes; ++i) {
        fprintf(cfo, "n %d %lf\n", i + NODE_IND_OFFSET, supply(i));
    }
    fclose(cfo);
    return;
}

void Graph::exportDatFile(const char *datFileOut) const
{
    FILE *dfo = fopen(datFileOut, "w");
    if (dfo == NULL) {
        printf("Unable to open file %s\n", datFileOut);
        exit(-1);
    }
    fprintf(dfo, "data;\n");

    fprintf(dfo, "set Nodes := ");
    for (int i = 0; i < numNodes; ++i) {
        fprintf(dfo, "%d,", i + NODE_IND_OFFSET);
    }
    fprintf(dfo, ";\n");

    fprintf(dfo, "set Arcs := ");
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNodes; ++j) {
            if (capacity(i, j) <= 0.0) { // Arc doesn't exist
                continue; 
            } // else arc exists
            fprintf(dfo, "(%d,%d),", i + NODE_IND_OFFSET, j + NODE_IND_OFFSET);
        }
    }
    fprintf(dfo, ";\n");

    fprintf(dfo, "param p := %d;\n", numEqualFlowSets);
    for (int r = 0; r < numEqualFlowSets; ++r) {
        fprintf(dfo, "set R[%d] := ", r + EQF_IND_OFFSET);
        for (int i = 0; i < numNodes; ++i) {
            for (int j = 0; j < numNodes; ++j) {
                if (capacity(i, j) <= 0.0) { // Arc doesn't exist
                    continue;
                } // else arc exists
                if (equalFlowIndex(i, j) == r) {
                    fprintf(dfo, "(%d,%d),", i + NODE_IND_OFFSET,
                                             j + NODE_IND_OFFSET);
                }
            }
        }
        fprintf(dfo, ";\n");
    }

    fprintf(dfo, "param: b := \n");
    for (int i = 0; i < numNodes; ++i) {
        fprintf(dfo, "%d %lf\n", i + NODE_IND_OFFSET, supply(i));
    }
    fprintf(dfo, ";\n");

    fprintf(dfo, "param u := \n");
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNodes; ++j) {
            if (capacity(i, j) <= 0.0) { // Arc doesn't exist
                continue;
            } // else arc exists
            fprintf(dfo, "%d %d %lf\n", i + NODE_IND_OFFSET, 
                    j + NODE_IND_OFFSET, capacity(i, j));
        }
    }
    fprintf(dfo, ";\n");

    fprintf(dfo, "param c := \n");
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNodes; ++j) {
            if (capacity(i, j) <= 0.0) { // Arc doesn't exist
                continue;
            } // else arc exists
            fprintf(dfo, "%d %d %lf\n", i + NODE_IND_OFFSET, 
                    j + NODE_IND_OFFSET, cost(i, j));
        }
    }
    fprintf(dfo, ";\n");

    fprintf(dfo, "param mu := \n");
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNodes; ++j) {
            if (capacity(i, j) <= 0.0) { // Arc doesn't exist
                continue;
            } // else arc exists
            fprintf(dfo, "%d %d %lf\n", i + NODE_IND_OFFSET, 
                    j + NODE_IND_OFFSET, multiplier(i, j));
        }
    }
    fprintf(dfo, ";\n");

    fprintf(dfo, "end;\n");
    fclose(dfo);
    return;
}

/*****************************************************************************/
/* Protected functions                                                       */
/*****************************************************************************/
void Graph::initializeDataStructures()
{
    // Create all nodes (uninitialized values, though)
    for (int i = 0; i < numNodes; ++i) {
        nodes.push_back(new GraphNode());
    }
    // Initialize matrices for storing arc information
    arcCosts.initialize(numNodes, 0.0);
    arcCapacities.initialize(numNodes, 0.0);
    arcMultipliers.initialize(numNodes, 0.0);
    eqFlowSetIndices.initialize(numNodes, -1);
    eqFlowNodeValues.initialize(numNodes, numEqualFlowSets + 1, 0.0);
    return;
}

// Add in self-loops, and place them in one big equal flow set
void Graph::addSelfLoops()
{
	for (int i = 1; i < numNodes; ++i)
	{
		// This is a little weird -- we have a negative multiplier on some
		// of the self-loops.  I don't know how the algorithm is going to 
		// react to that
		if (supply(i) != 0)
		{
			arcCosts.set(i, i, bigM);
			arcCapacities.set(i, i, std::numeric_limits<double>::max());
			arcMultipliers.set(i, i, -supply(i) + 1);
			eqFlowSetIndices.set(i, i, numEqualFlowSets);
			eqFlowNodeValues.increment(i, numEqualFlowSets, supply(i));
		}
	}
}

void Graph::readColFile(const char *filename)
{
    std::ifstream inFile;
    string nextLine;

    inFile.open(filename);
    if (!inFile) {
        printf("Unable to open file: %s\n", filename);
        exit(1);
    }   

    // Read from file
    int arcsRead = 0;
    while (!inFile.eof()) {
        fflush(stdout);    
        getline(inFile, nextLine);
        if (nextLine.empty()) {
            continue;
        }
        switch (nextLine[0]) {
        case 'p':
            processProblemLine(nextLine);
            break;
        case 'e':
            processEdgeLine(nextLine);
            break;
        case 'a':
            ++arcsRead;
            processArcLine(nextLine);
            break;
        case 'n':
            processNodeLine(nextLine);
            break;
        case 'c':
        default:
            // Unrecognized lines treated as comments
            processCommentLine(nextLine);
            break;
        }
    }
    // Close file
    inFile.close();

    if (arcsRead != numArcs) {
        printf("Error: Number of arcs specified does not match the number of "
               "arcs in the file\n");
        exit(1);
    }

    return;
}

void Graph::readGSFile(const char *filename)
{
    std::ifstream inFile;
    string nextLine;

    // First pass to count number of node-arc lines and equal flow sets
    inFile.open(filename);
    if (!inFile) {
        printf("Unable to open file: %s\n", filename);
        exit(1);
    }
    while (!inFile.eof()) {
        fflush(stdout);    
        getline(inFile, nextLine);
        if ((nextLine.empty()) || (nextLine[0] != '(')) {
            continue;
        }
        int found = nextLine.find(":");
        if (found != string::npos) { // Arc
            ++numNodes;
        } else { // Equal flow set
            ++numEqualFlowSets;
        }
    }
    inFile.close();

    // Now intialize graph data structures
    initializeDataStructures();

    // Second pass to actually process the file
    int eqFlowNum = 0;
    inFile.open(filename);
    if (!inFile) {
        printf("Unable to open file: %s\n", filename);
        exit(1);
    }
    while (!inFile.eof()) {
        fflush(stdout);    
        getline(inFile, nextLine);
        if ((nextLine.empty()) || (nextLine[0] != '(')) {
            continue;
        }
        int found = nextLine.find(":");
        if (found != string::npos) { // Arc
            parseGSArcLine(nextLine);
        } else { // Equal flow set
            parseGSEqualFlowLine(nextLine, eqFlowNum);
            ++eqFlowNum;
        }
    }
    inFile.close();

    return;
}

/*****************************************************************************/
/* Functions for handling .col files                                         */
/*****************************************************************************/
void Graph::processProblemLine(string nextLine)
{
    int result;
    char tempString[10];

//    printf("Problem line: %s\n", nextLine.c_str());
    result = sscanf(nextLine.c_str(), "p %s %d %d %d \n", tempString, 
                    &numNodes, &numArcs, &numEqualFlowSets);
    if (result != 4) {
        printf("Improperly formatted problem line\n");
        exit(1);
    }
    initializeDataStructures();
    return;
}

void Graph::processEdgeLine(string nextLine)
{
//    printf("Edge line: %s\n", nextLine.c_str());
//    printf("Edges are only allowed in undirected graphs; ignoring...\n");
    processArcLine(nextLine);
    return;
}

void Graph::processArcLine(string nextLine)
{
    int result;
    int n1, n2, eqFlowInd;
    double tempLB, tempUB, tempCost, tempMult;

//    printf("Arc line: %s\n", nextLine.c_str());
    result = sscanf(nextLine.c_str(), "a %d %d %lf %lf %lf %lf %d \n", 
                &n1, &n2, &tempLB, &tempUB, &tempCost, &tempMult, &eqFlowInd);
    if (result != 7) {
        printf("Improperly formatted arc line\n");
        exit(1);
    }

    if (arcCapacities.get(n1-1, n2-1) > 0.0) {
        printf("Warning! Duplicate arc detected, behaviour is undefined.\n");
    }
//    if (arcCapacities.get(n2-1, n1-1) > 0.0) {
//        printf("Reverse arc detected\n");
//    }
//    if (n1 == n2) {
//        printf("Self-loop detected in input graph.\n");
//    }

    arcCosts.set(n1 - NODE_IND_OFFSET, n2 - NODE_IND_OFFSET, tempCost);
    arcCapacities.set(n1 - NODE_IND_OFFSET, n2 - NODE_IND_OFFSET, tempUB);
    arcMultipliers.set(n1 - NODE_IND_OFFSET, n2 - NODE_IND_OFFSET, tempMult);
    eqFlowSetIndices.set(n1 - NODE_IND_OFFSET, n2 - NODE_IND_OFFSET, 
                         eqFlowInd - EQF_IND_OFFSET);

    if (eqFlowInd - EQF_IND_OFFSET >= 0) {
        eqFlowNodeValues.increment(n1 - NODE_IND_OFFSET, 
                                   eqFlowInd - EQF_IND_OFFSET,  
                                    1.0);
        eqFlowNodeValues.increment(n2 - NODE_IND_OFFSET, 
                                   eqFlowInd - EQF_IND_OFFSET, 
                                   -1.0 * tempMult);
    }

    if (maxArcCost < tempCost) {
        maxArcCost = tempCost;
    }
    if (maxArcCapacity < tempUB) {
        maxArcCapacity = tempUB;
    }
    return;
}

void Graph::processNodeLine(string nextLine)
{
    int result;
    int node1;
    double tempSupply;

//    printf("Node line: %s\n", nextLine.c_str());
    result = sscanf(nextLine.c_str(), "n %d %lf \n", &node1, &tempSupply);
    if (result != 2) {
        printf("Improperly formatted node line\n");
        exit(1);
    }
    nodes[node1 - NODE_IND_OFFSET]->id = node1;
    nodes[node1 - NODE_IND_OFFSET]->index = node1 - NODE_IND_OFFSET;
    nodes[node1 - NODE_IND_OFFSET]->supply = tempSupply;
    totalSupply += tempSupply; // Not really sure why we need this
    return;
}

void Graph::processCommentLine(string nextLine)
{
//    printf("Comment line: %s\n", nextLine.c_str());
    return;
}

/*****************************************************************************/
/* Functions for handling .gs files                                          */
/*****************************************************************************/
void Graph::parseGSArcLine(string nextLine)
{
    int result, beginInd, endInd;
    int node1, node2;
    double tempSupply, tempCost, tempUB, tempMult; 

//    printf("Node-Arc line: %s\n", nextLine.c_str());

    // Identify node
    beginInd = nextLine.find("(");
    endInd = nextLine.find(")");
    string nodePortion = nextLine.substr(beginInd, endInd - beginInd + 1);
    nextLine = nextLine.substr(endInd + 1);
//    printf("Identified the node %s\n", nodePortion.c_str());
//    printf("And the remainder is %s\n", nextLine.c_str());
    result = sscanf(nodePortion.c_str(), "(%d %lf)", &node1, &tempSupply);
    nodes[node1 - NODE_IND_OFFSET]->id = node1;
    nodes[node1 - NODE_IND_OFFSET]->index = node1 - NODE_IND_OFFSET;
    nodes[node1 - NODE_IND_OFFSET]->supply = tempSupply;
    totalSupply += tempSupply;

    while (!(nextLine.empty())) {
        beginInd = nextLine.find("(");
        endInd = nextLine.find(")");
        if (beginInd == string::npos) {
            break;
        }
        string arcPortion = nextLine.substr(beginInd, endInd - beginInd + 1);
        nextLine = nextLine.substr(endInd + 1);
//        printf("Identified the arc %s\n", arcPortion.c_str());
//        printf("And the remainder is %s\n", nextLine.c_str());
        result = sscanf(arcPortion.c_str(), "(%d %lf %lf %lf)", &node2, 
            &tempUB, &tempCost, &tempMult);

        if (capacity(node1 - NODE_IND_OFFSET, node2 - NODE_IND_OFFSET) <= 0.0) {
            ++numArcs;
        } else {
            printf("Duplicate arc detected\n");
        }

        arcCosts.set(node1 - NODE_IND_OFFSET, node2 - NODE_IND_OFFSET, 
                        tempCost);
        arcCapacities.set(node1 - NODE_IND_OFFSET, node2 - NODE_IND_OFFSET, 
                        tempUB);
        arcMultipliers.set(node1 - NODE_IND_OFFSET, node2 - NODE_IND_OFFSET, 
                        tempMult);
        eqFlowSetIndices.set(node1 - NODE_IND_OFFSET, node2 - NODE_IND_OFFSET, 
                             -1);

        if (maxArcCost < tempCost) {
            maxArcCost = tempCost;
        }
        if (maxArcCapacity < tempUB) {
            maxArcCapacity = tempUB;
        }
    }

    return;
}

void Graph::parseGSEqualFlowLine(string nextLine, int eqFlowNum)
{
    int result, beginInd, endInd;
    int node1, node2;

//    printf("EqFlowSet line: %s\n", nextLine.c_str());

    while (!(nextLine.empty())) {
        beginInd = nextLine.find("(");
        endInd = nextLine.find(")");
        if (beginInd == string::npos) {
            break;
        }
        string arcPortion = nextLine.substr(beginInd, endInd - beginInd + 1);
        nextLine = nextLine.substr(endInd + 1);
//        printf("Identified the arc %s\n", arcPortion.c_str());
//        printf("And the remainder is %s\n", nextLine.c_str());
        result = sscanf(arcPortion.c_str(), "(%d %d)", &node1, &node2);

        eqFlowSetIndices.set(node1 - NODE_IND_OFFSET, node2 - NODE_IND_OFFSET, 
                             eqFlowNum);
        eqFlowNodeValues.increment(node1 - NODE_IND_OFFSET, eqFlowNum, 1.0);
        eqFlowNodeValues.increment(node2 - NODE_IND_OFFSET, eqFlowNum, -1.0 * 
            multiplier(node1 - NODE_IND_OFFSET, node2 - NODE_IND_OFFSET));
    }

    return; 
}

