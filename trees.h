/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-06-16                                                          */
/* File: trees.h                                                             */
/* Description:                                                              */
/*   Contains design details for the tree structures used in the GenNetEq    */
/*   algorithm.                                                              */
/*****************************************************************************/
#ifndef TREES_H
#define TREES_H

// Required include's
#include <vector>
using std::vector;

// Forward declarations
class Node;
class ArcVar;

class Tree;
typedef Tree TreeTypeI;
typedef Tree TreeTypeII;

class Tree 
{
  public:
  // Public static class methods
    static void setStaticVariables(vector<Node> *gneNodes, 
                                   vector<TreeTypeI *> *gneSetBtI, 
                                   vector<TreeTypeII *> *gneSetBtII);

    static void printTrees(bool detailed = false);
    static void printTree(Tree *tX, bool detailed = false);
    static void printTreeNodes(bool detailed = false);
    static void printTreeNode(int nodeID, bool detailed = false);

    // Used for debugging trees
    static void checkTrees();

    // Used when inserting or removing arcs
    static void insertIntoTrees(ArcVar *av);
    static void checkForStrandedNodes();
    static void removeFromTrees(ArcVar *av);
    static void updateTrees();

  // INLINED FUNCTIONS
    // Miscellaneous getters
    int getTreeID();
    tree_type getTreeType();
    int getRoot();
    ArcVar *getExtraArc();
    vector<int> *getNodesInTree();
    vector<ArcVar *> *getArcsInTree();

	int getSize() { return nodesInTree.size(); }

    // Getters for thread indices
    void resetRootToLeafTraversal();
    void resetLeafToRootTraversal();
    int getNextThread();
    int getPrevThread();

    // Dynamic tree constructor and destructor
    ~Tree();

  protected:
    Tree();

  // Protected dynamic class variables
    int treeID;
    tree_type treeType;
    int root;
    bool upToDate;
    ArcVar *extraArc;
    vector<int> nodesInTree;
    vector<ArcVar *> arcsInTree;

    // Needed for computing flows and potentials using thread indices
    vector<int> threads; 
    int threadPos;

  // Protected static class variables
    static int numNodes;
    static vector<Node> *nodes;
    static vector<TreeTypeI *> *setBtI;
    static vector<TreeTypeII *> *setBtII;

  // Protected static class methods
    // Miscellaneous
    static Tree *createNewTree(tree_type tType, int rNode = -1);    // inline
    static void initializeTree(Tree *tX, int tID, tree_type tType); // inline
    static void setTreeIDAndType(Tree *tX, int newID, tree_type newType);
    static void removeTree(Tree *tX);                               // inline
    static void convertTItoTII(TreeTypeI *tI);                      // inline
    static void convertTIItoTI(TreeTypeII *tII, ArcVar *av);        // inline

    // Used when inserting arcs
    static void mergeTxAndTII(Tree *tX, TreeTypeII *tII, ArcVar *av);
    static void expandTreeUsingArc(Tree *tX, ArcVar *av);           // inline

    // Used when removing arcs
    static void splitTree(Tree *tX, ArcVar *removedArc);

  // Protected dynamic tree methods for updating indices
    void initializeTreeIndices();

  private:
    // Nothing
};

/*****************************************************************************/
/* Dynamic tree class inline functions                                       */
/*****************************************************************************/
inline
int Tree::getTreeID()
{
    return treeID;
}

inline
tree_type Tree::getTreeType()
{
    return treeType;
}

inline
int Tree::getRoot()
{
    return root;
}

inline
ArcVar *Tree::getExtraArc()
{
    return extraArc;
}

inline
vector<int> *Tree::getNodesInTree()
{
    return &nodesInTree;
}

inline
vector<ArcVar *> *Tree::getArcsInTree()
{
    return &arcsInTree;
}

// Getters for computing flows and potentials using thread indices
inline
void Tree::resetRootToLeafTraversal()
{
    threadPos = 0;
    return;
}

inline
void Tree::resetLeafToRootTraversal()
{
    threadPos = threads.size() - 1;
    return;
}

inline 
int Tree::getNextThread()
{
    return (threadPos == threads.size()) ? threads[0] : threads[threadPos++];
}

inline
int Tree::getPrevThread()
{
    return (threadPos < 0) ? threads[0] : threads[threadPos--];
}

#endif // TREES_H

