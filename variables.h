/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-06-16                                                          */
/* File: variables.h                                                         */
/* Description:                                                              */
/*   Contains design details for network simplex algorithm variables.        */
/*****************************************************************************/
#ifndef VARIABLES_H
#define VARIABLES_H

// Required include's
#include <cstdio>
#include <cstdlib>

#include <vector>
using std::vector;

// Global Variable variables
const int THETA = -1;
const int NULL_THETA = -2;
const int NULL_TREE_ID = -1;

enum var_loc {B_var, L_var, U_var, NIL};
enum var_type { ARC_VAR, EQF_VAR };
enum tree_type { NULL_TREE_TYPE, TYPE_I, TYPE_II };

// Forward Declarations
class ArcVar;

struct ThetaEquation
{
    double a0;
    double a1;
    int thetaID;
};

class Node
{
  public:
    Node()
    {
        potentialWithEq.a0 = 0.0;
        potentialWithEq.a1 = 0.0;
        potentialWithEq.thetaID = NULL_THETA;

        supplyWithEq.a0 = 0.0;
        supplyWithEq.a1 = 0.0;
        supplyWithEq.thetaID = NULL_THETA;

        memberOfTreeID = NULL_TREE_ID;
        memberOfTreeType = NULL_TREE_TYPE;

        depth = -1;
        thread = -1;
        pred = -1;
        predArc = NULL;
    }
    ~Node() { /* do nothing */ };

    int id; 

    double potential;
    double supply;

    ThetaEquation potentialWithEq;
    ThetaEquation supplyWithEq;
    vector<double> supplyWithThetas;

    // Member variables
    int memberOfTreeID;
    tree_type memberOfTreeType;

    vector<ArcVar *> incidentArcs;

    // NOTE: Technically we don't need depth and thread, but we'll maintain 
    // them anyways
    int depth;
    int thread;
    int pred;
    ArcVar *predArc;

  protected:
    // Nothing

  private:
    // Nothing
};

class Variable
{
  public:
    Variable()
    {
        setLoc = NIL;
        setInd = -1;
        flow = 0.0;
        deltaFlow = 0.0;
        redCost = 0.0;

        flowWithEq.a0 = 0.0;
        flowWithEq.a1 = 0.0;
        flowWithEq.thetaID = NULL_THETA;
    }
    virtual ~Variable() { /* do nothing */ };

    int varIndex;
    double LB;
    double UB;
    double cost;
    double actualCost;
    enum var_type varType;
    bool isArtificial;

    // Variable location information
    var_loc setLoc;
    int setInd;

    // Flow stuff
    double flow;
    double deltaFlow;
    double redCost;

    ThetaEquation flowWithEq;
    vector<double> flowWithThetas;

    virtual void print(bool detailed = false) = 0;

  protected:
    // Nothing

  private:
    // Nothing
};

class ArcVar : public Variable
{
  public:
    ArcVar() { /* do nothing */ };
    virtual ~ArcVar() { /* do nothing */ };

    int tail;
    int head;
    double mult;

    virtual void print(bool detailed = false)
    {
        printf("%d: (%d,%d)", varIndex, tail, head);
        if (detailed) {
            if (isArtificial) {
                printf("-F");
            }
            printf(" [%f,%f] %f, %f", LB, UB, cost, mult);
        }
        printf("\n");
        return;
    }

  protected:
    // Nothing

  private:
    // Nothing
};

class EqFlowVar : public Variable
{
  public:
    EqFlowVar() { /* do nothing */ };
    virtual ~EqFlowVar() { /* do nothing */ };

    int eqFlowIndex;
    int numOrigArcs;

    vector<double> nodeVals;

    virtual void print(bool detailed = false)
    {
        printf("%d: R_{%d}", varIndex, eqFlowIndex);
        if (detailed) {
            printf(" [%f,%f] %f, %d", LB, UB, cost, numOrigArcs);
        }
        printf("\n");
        return;
    }

  protected:
    // Nothing

  private:
    // Nothing
};

// Variable Comparator
struct varcomp 
{
    bool operator() (Variable *v1, Variable *v2) 
    {
        return (v1->varIndex < v2->varIndex);
    }
};

struct PivotVars 
{
	int eVarInd;
	int lVarInd;
};

struct VarWithViol
{
    Variable *var;
    double violation;
};

// VarWithViol Comparator
struct violcomp 
{
    bool operator() (VarWithViol *v1, VarWithViol *v2) 
    {
        return (v1->violation < v2->violation);
    }
};

#endif // VARIABLES_H

