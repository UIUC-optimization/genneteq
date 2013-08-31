/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-06-16                                                          */
/* File: graph.h                                                             */
/* Description:                                                              */
/*   Contains design details for a custom graph class.                       */
/*****************************************************************************/
#ifndef GRAPH_H
#define GRAPH_H

// Required include's
#include <vector>
using std::vector;
#include <string>
using std::string;

const int NODE_IND_OFFSET = 1;
const int EQF_IND_OFFSET = 1;

// Forward Declarations
struct Config;
template <typename T> class Matrix;

struct GraphNode {
    int id;
    int index;
    double supply;
};

class Graph {
  public:
    Graph();
    ~Graph();

    void initialize(Config *conf);
    void print() const;
    void exportColFile(const char *colFileOut) const;
    void exportDatFile(const char *datFileOut) const;

    // Inlined functions
    int getNumNodes() const;
    int getNumArcs() const;
    int getNumEqualFlowSets() const;
    double getBigM() const;

    double supply(int i) const;
    double cost(int i, int j) const;
    double capacity(int i, int j) const;
    double multiplier(int i, int j) const;
    int equalFlowIndex(int i, int j) const;
    double equalFlowNodeValue(int i, int r) const;

  protected:
    // Variables
    int numNodes;
    int numArcs;
    int numEqualFlowSets;

    double totalSupply;
    double maxArcCost;
    double maxArcCapacity;
    double bigM;

    vector<GraphNode *> nodes;
    Matrix<double> arcCosts;
    Matrix<double> arcCapacities;
    Matrix<double> arcMultipliers;
    Matrix<int> eqFlowSetIndices;
    Matrix<double> eqFlowNodeValues;

    // Functions
    void initializeDataStructures();
	void addSelfLoops();
	void addSelfLoop(int i);

    void readColFile(const char *filename);
    void readGSFile(const char *filename);

    void processProblemLine(string nextLine);
    void processEdgeLine(string nextLine);
    void processArcLine(string nextLine);
    void processNodeLine(string nextLine);
    void processCommentLine(string nextLine);

    void parseGSArcLine(string nextLine);
    void parseGSEqualFlowLine(string nextLine, int eqFlowNum);


  private:
    // Nothing
};

/*****************************************************************************/
/* Inline function declarations                                              */
/*****************************************************************************/
inline
int Graph::getNumNodes() const 
{
    return numNodes;
}

inline
int Graph::getNumArcs() const
{
    return numArcs;
}

inline
int Graph::getNumEqualFlowSets() const
{
    return numEqualFlowSets;
}

inline
double Graph::getBigM() const
{
    return bigM;
}

inline
double Graph::supply(int i) const 
{
    return nodes[i]->supply;
}

inline
double Graph::cost(int i, int j) const 
{
    return arcCosts.get(i, j);
}

inline
double Graph::capacity(int i, int j) const 
{
    return arcCapacities.get(i, j);
}

inline
double Graph::multiplier(int i, int j) const
{
    return arcMultipliers.get(i, j);
}

inline
int Graph::equalFlowIndex(int i, int j) const
{
    return eqFlowSetIndices.get(i, j);
}

inline
double Graph::equalFlowNodeValue(int i, int r) const
{
    return eqFlowNodeValues.get(i, r);
}

#endif // GRAPH_H

