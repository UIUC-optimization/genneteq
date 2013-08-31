/*****************************************************************************/
/* Author: Jason Sauppe                                                      */
/* Date: 2010-06-16                                                          */
/* File: trees.cpp                                                           */
/* Description:                                                              */
/*   Contains implementation details for the tree structures used in the     */
/*   GenNetEq algorithm.                                                     */
/*****************************************************************************/
#include "main.h"
#include "variables.h"
#include "trees.h"

#include <cstdio>
#include <cstdlib>

#include <vector>
using std::vector;
#include <stack>
using std::stack;

// Declaration of static variables
int Tree::numNodes;
vector<Node> *Tree::nodes;
vector<TreeTypeI *> *Tree::setBtI;
vector<TreeTypeII *> *Tree::setBtII;

/*****************************************************************************/
/* Dynamic tree constructor and destructor                                   */
/*****************************************************************************/
Tree::Tree()
{
    // Default constructor, do nothing
}

Tree::~Tree()
{
    // Clean-up, nothing to clean 
}

/*****************************************************************************/
/* Miscellaneous static tree class functions                                 */
/*****************************************************************************/
// Running time: O(1)
void Tree::setStaticVariables(vector<Node> *gneNodes, 
                              vector<TreeTypeI *> *gneSetBtI,
                              vector<TreeTypeII *> *gneSetBtII)
{
    Tree::nodes = gneNodes;
    Tree::setBtI = gneSetBtI;
    Tree::setBtII = gneSetBtII;
    Tree::numNodes = Tree::nodes->size();
    return;
}

// Running time: O(n)
void Tree::printTrees(bool detailed)
{
    printf("Trees are:\n");
    for (int tIind = 0; tIind < (*setBtI).size(); ++tIind) {
        printTree((*setBtI)[tIind], detailed);
    }
    for (int tIIind = 0; tIIind < (*setBtII).size(); ++tIIind) {
        printTree((*setBtII)[tIIind], detailed);
    }
    return;
}

// Running time: O(|Ti|)
void Tree::printTree(Tree *tX, bool detailed)
{
    printf("Tree %d of type %d contains %d nodes and %d arcs (root = %d) %d\n", 
            tX->treeID, tX->treeType, tX->nodesInTree.size(), 
            tX->arcsInTree.size(), tX->root, tX->upToDate);
    if (detailed) {
        printf("  Nodes: [");
        for (int nodeInd = 0; nodeInd < tX->nodesInTree.size(); ++nodeInd) {
            printf("%d,", (tX->nodesInTree)[nodeInd]);
        }
        printf("]\n  Arcs: [");
        for (int arcInd = 0; arcInd < tX->arcsInTree.size(); ++arcInd) {
            printf("(%d,%d)", (tX->arcsInTree)[arcInd]->tail, 
                              (tX->arcsInTree)[arcInd]->head);
            if ((tX->arcsInTree)[arcInd]->isArtificial) {
                printf("-F");
            }
            printf(",");
        }
        printf("]\n");
        if (tX->treeType == TYPE_I) {
            printf("  Ex.Arc: (%d,%d)", tX->extraArc->tail, tX->extraArc->head);
            if (tX->extraArc->isArtificial) {
                printf("-F");
            }
            printf("\n");
        }
    }
    return;
}

// Running time: O(n)
void Tree::printTreeNodes(bool detailed)
{
    for (int i = 0; i < numNodes; ++i) {
        printTreeNode(i, detailed);
    }
    return;
}

// Running time: O(d_Ti(u))
void Tree::printTreeNode(int nodeID, bool detailed)
{
    Node &tn = (*Tree::nodes)[nodeID];
    printf("  %4d: Tree:=%d,%d; D:=%d; T:=%d; P:=%d\n", 
            nodeID, tn.memberOfTreeID, tn.memberOfTreeType, 
            tn.depth, tn.thread, tn.pred);
    if (detailed) {
        if (tn.predArc == NULL) {
            printf("    PredArc:= NULL\n");
        } else {
            printf("    PredArc:= (%d,%d)\n", tn.predArc->tail, 
                                              tn.predArc->head);
        }
        printf("    Incident Arcs: [");
        for (int arcInd = 0; arcInd < tn.incidentArcs.size(); ++arcInd) {
            printf("(%d,%d)", tn.incidentArcs[arcInd]->tail, 
                              tn.incidentArcs[arcInd]->head);
            if (tn.incidentArcs[arcInd]->isArtificial) {
                printf("-F");
            }
            printf(",");
        }
        printf("]\n");
    }
    return;
}

/*****************************************************************************/
/* Static tree debugging functions                                           */
/*****************************************************************************/
// Running time: O(# trees)
void Tree::checkTrees()
{
    for (int tIind = 0; tIind < (*setBtI).size(); ++tIind) {
        TreeTypeI *tI = (*setBtI)[tIind];
        if (tI->treeID != tIind) {
            printf("Error: Tree ID is not consistent with index\n");
            exit(-1);
        }
        if (tI->extraArc == NULL) {
            printf("Error: Type I tree is missing extra arc\n");
            exit(-1);
        }
        if (tI->nodesInTree.size() != tI->arcsInTree.size() + 1) {
            printf("Error: Nodes (%d) != Arcs + 1 (%d)\n", 
                tI->nodesInTree.size(), tI->arcsInTree.size());
            exit(-1);
        }
    }
    for (int tIIind = 0; tIIind < (*setBtII).size(); ++tIIind) {
        TreeTypeII *tII = (*setBtII)[tIIind];
        if (tII->treeID != tIIind) {
            printf("Error: Tree ID is not consistent with index\n");
            exit(-1);
        }
        if (tII->extraArc != NULL) {
            printf("Error: Type II tree has extra arc\n");
            exit(-1);
        }
        if (tII->nodesInTree.size() != tII->arcsInTree.size() + 1) {
            printf("Error: Nodes (%d) != Arcs + 1 (%d)\n", 
                tII->nodesInTree.size(), tII->arcsInTree.size());
            exit(-1);
        }
    }
    return;
}

/*****************************************************************************/
/* Static tree constructors and destructors                                  */
/*****************************************************************************/
// Running time: O(1)
inline
Tree *Tree::createNewTree(tree_type tType, int rNode)
{
    Tree *t = new Tree();
    if (tType == TYPE_I) {
        initializeTree(t, (*Tree::setBtI).size(), TYPE_I);
        (*Tree::setBtI).push_back(t);
    } else { // tType == TYPE_II
        initializeTree(t, (*Tree::setBtII).size(), TYPE_II);
        (*Tree::setBtII).push_back(t);
    }

    // Inserting node into t and setting it as root
    if (rNode >= 0) {
        t->root = rNode;
        t->nodesInTree.push_back(rNode);
        (*Tree::nodes)[rNode].memberOfTreeID = t->treeID;
        (*Tree::nodes)[rNode].memberOfTreeType = t->treeType;
    }

    return t;
}

// Running time: O(1)
inline
void Tree::initializeTree(Tree *tX, int tID, tree_type tType)
{
    tX->treeID = tID;
    tX->treeType = tType;
    tX->root = -1;
    tX->upToDate = false;
    tX->extraArc = NULL;
    return;
}

// Running time: O(|Ti|)
void Tree::setTreeIDAndType(Tree *tX, int newID, tree_type newType)
{
    tX->treeID = newID;
    tX->treeType = newType;
    for (int tInd = 0; tInd < tX->nodesInTree.size(); ++tInd) {
        int nodeID = (tX->nodesInTree)[tInd];
        (*Tree::nodes)[nodeID].memberOfTreeID = newID;
        (*Tree::nodes)[nodeID].memberOfTreeType = newType;
    }
    return;
}

// Running time: O(|Tj|)
inline
void Tree::removeTree(Tree *tX)
{
    // To remove a tree, simply copy the last tree of appropriate type into 
    // the spot occupied by this tree and then erase the last tree in the 
    // appropriate vector
    int oldIDtX = tX->treeID;
    if (tX->treeType == TYPE_I) {
        if (oldIDtX < (*Tree::setBtI).size() - 1) {
            TreeTypeI *lastTI = (*Tree::setBtI).back();
            setTreeIDAndType(lastTI, oldIDtX, TYPE_I);
            (*Tree::setBtI)[oldIDtX] = lastTI;
        }
        (*Tree::setBtI).pop_back();
    } else { // tX->treeType == TYPE_II
        if (oldIDtX < (*Tree::setBtII).size() - 1) {
            TreeTypeII *lastTII = (*Tree::setBtII).back();
            setTreeIDAndType(lastTII, oldIDtX, TYPE_II);
            (*Tree::setBtII)[oldIDtX] = lastTII;
        }
        (*Tree::setBtII).pop_back();
    }
    return;
}

// Running time: O(|Ti| + |Tj|) = O(n)
inline
void Tree::convertTItoTII(TreeTypeI *tI)
{
    int oldIDtI = (*tI).getTreeID();
    int newIDtII = (*Tree::setBtII).size();
    // Remove the old type II tree [O(|Tj|)]
    Tree::removeTree(tI);
    // Convert old type I tree to a type II tree without extra arc [O(|Ti|)]
    setTreeIDAndType(tI, newIDtII, TYPE_II);
    tI->extraArc = NULL;
    (*Tree::setBtII).push_back(tI);
    return;
}

// Running time: O(|Ti| + |Tj|) = O(n)
inline
void Tree::convertTIItoTI(TreeTypeII *tII, ArcVar *av)
{
    int oldIDtII = (*tII).getTreeID();
    int newIDtI = (*Tree::setBtI).size();
    // Remove the old type II tree [O(|Tj|)]
    Tree::removeTree(tII);
    // Convert old type II tree to a type I tree with extra arc [O(|Ti|)]
    setTreeIDAndType(tII, newIDtI, TYPE_I);
    tII->extraArc = av;
    (*Tree::setBtI).push_back(tII);
    return;
}

/*****************************************************************************/
/* Insert Arc Methods                                                        */
/*****************************************************************************/
// Running time: O(n)
void Tree::insertIntoTrees(ArcVar *av)
{
    int u = av->tail;
    int v = av->head;
    int tUID = (*Tree::nodes)[u].memberOfTreeID;
    tree_type tUType = (*Tree::nodes)[u].memberOfTreeType;
    int tVID = (*Tree::nodes)[v].memberOfTreeID;
    tree_type tVType = (*Tree::nodes)[v].memberOfTreeType;

    // DEBUG:
//    printf("Arc being inserted is: (%d,%d)\n", u, v);
//    fflush(stdout);

    if ((tUID != NULL_TREE_ID) && (tVID != NULL_TREE_ID)) {
        // Both nodes in arc belong to trees, so need to merge or create cycle
//        printf("  Working with existing trees\n");
        if ((tUID == tVID) && (tUType == tVType)) {
            // Both nodes belong to the same tree, so a cycle is created. 
            // The tree better be a type II tree, otherwise it will have two 
            // cycles which is not allowed. [O(n)]
            if (tUType == TYPE_I) {
                printf("Error: Insertion adds cycle to type I tree\n");
                exit(-1);
            }
            Tree *tUV = (*Tree::setBtII)[tUID];
            Tree::convertTIItoTI(tUV, av);
        } else { 
            // Nodes belong to two different trees so need to merge them [O(n)]
            if ((tUType == TYPE_II) &&  (tVType == TYPE_II)) {
                mergeTxAndTII((*Tree::setBtII)[tUID], 
                              (*Tree::setBtII)[tVID], av);
            } else if ((tUType == TYPE_I) && (tVType == TYPE_II)) {
                mergeTxAndTII((*Tree::setBtI)[tUID], 
                              (*Tree::setBtII)[tVID], av);
            } else if ((tUType == TYPE_II) && (tVType == TYPE_I)) {
                mergeTxAndTII((*Tree::setBtI)[tVID], 
                              (*Tree::setBtII)[tUID], av);
            } else { // ((tUType == TYPE_I) && (tVType == TYPE_I))
                printf("Error: Insertion merges two type I trees\n");
                exit(-1);
            }
        }
    } else { // Grow tree the old-fashioned way [O(1)]
        // At least one of the nodes does not yet belong in a tree, so we need 
        // to either append the arc to an existing tree or create a new tree.
        // Note that this *cannot* create any cycles except for self-loops
//        printf("  Growing a tree the old-fashioned way\n");
        Tree *tUV;
        if ((tUID == NULL_TREE_ID) && (tVID == NULL_TREE_ID)) {
            // Need to create a new tree for arc (u,v) [O(1)]
            if (u == v) {
                tUV = Tree::createNewTree(TYPE_I, u);
            } else {
                tUV = Tree::createNewTree(TYPE_II, u);
            }
        } else if ((tUID != NULL_TREE_ID) && (tVID == NULL_TREE_ID)) {
            if (tUType == TYPE_I) {
                tUV = (*Tree::setBtI)[tUID];
            } else { // tUType == TYPE_II
                tUV = (*Tree::setBtII)[tUID];
            }
        } else { // ((tUID == NULL_TREE_ID) && (tVID != NULL_TREE_ID)) {
            if (tVType == TYPE_I) {
                tUV = (*Tree::setBtI)[tVID];
            } else { // tVType == TYPE_II
                tUV = (*Tree::setBtII)[tVID];
            }
        } 
        expandTreeUsingArc(tUV, av); // [O(1)]
    }

    return;
}

// Running time: O(|Ti| + |Tj|) = O(n)
void Tree::mergeTxAndTII(Tree *tX, TreeTypeII *tII, ArcVar *av)
{
    // First remove the type II tree from the list of type II trees [O(|Tj|)]
    Tree::removeTree(tII);

    // Update tree node information for nodes in tII [O(|Ti|)]
    setTreeIDAndType(tII, tX->treeID, tX->treeType);

    // Copy nodes from tII into tX [O(|Ti|)]
    vector<int> *tIINodes = (*tII).getNodesInTree();
    for (int tIINodeInd = 0; tIINodeInd < (*tIINodes).size(); ++tIINodeInd) {
        tX->nodesInTree.push_back((*tIINodes)[tIINodeInd]);
    }
    // Copy arcs from tII into tX [O(|Ti|)]
    vector<ArcVar *> *tIIArcs = (*tII).getArcsInTree();
    for (int tIIArcInd = 0; tIIArcInd < (*tIIArcs).size(); ++tIIArcInd) {
        tX->arcsInTree.push_back((*tIIArcs)[tIIArcInd]);
    }

    // Add arc av to tX and update incident arcs for nodes u and v
    tX->arcsInTree.push_back(av);
    (*Tree::nodes)[av->tail].incidentArcs.push_back(av);
    (*Tree::nodes)[av->head].incidentArcs.push_back(av);

    tX->upToDate = false;

    // Now we want to delete tII, since we don't need it any more
    delete tII;

    return;
}

// Running time: O(1)
inline
void Tree::expandTreeUsingArc(Tree *tX, ArcVar *av)
{
    // NOTE: This expects that arcs are added to a tree in such a way so that 
    // the tree is always connected; i.e., each arc inserted must contain a 
    // node already in the tree, and the other node on the arc cannot belong 
    // to any other tree (except this one in the case of self-loops). If this 
    // is not the case, then the two trees should be merged when the arc is 
    // added. Note that this is handled by the insertIntoTrees routine. 
    int u = av->tail;
    int v = av->head;
    Node &tnU = (*Tree::nodes)[u];
    Node &tnV = (*Tree::nodes)[v];

    if (tnU.memberOfTreeID != tX->treeID) { // u not in tree, so v is
        tX->nodesInTree.push_back(u);
        tnU.memberOfTreeID = tX->treeID;
        tnU.memberOfTreeType = tX->treeType;
    } else if (tnV.memberOfTreeID != tX->treeID) { // v not in tree, so u is
        tX->nodesInTree.push_back(v);
        tnV.memberOfTreeID = tX->treeID;
        tnV.memberOfTreeType = tX->treeType;
    } // else both u and v are in tree, so this must be a self-loop

    // Now both nodes are stored in tree, so update incident arcs and adjacent 
    // tree nodes for them, assuming arc is not a self-loop
    if (u != v) {
        tX->arcsInTree.push_back(av);
        tnU.incidentArcs.push_back(av);
        tnV.incidentArcs.push_back(av);
    } else { 
        tX->extraArc = av;
    }

    // Any time this method is called, indices are destroyed
    tX->upToDate = false;

    return;
}

/*****************************************************************************/
/* Insert Node Methods                                                       */
/*****************************************************************************/
// Running time: O(n)
void Tree::checkForStrandedNodes()
{
    for (int i = 0; i < Tree::numNodes; ++i) {
        if ((*Tree::nodes)[i].memberOfTreeID == NULL_TREE_ID) {
            createNewTree(TYPE_II, i);
        }
    }
    return;
}

/*****************************************************************************/
/* Remove Arc Methods                                                        */
/*****************************************************************************/
// Running time: O(n)
void Tree::removeFromTrees(ArcVar *av)
{
    int u = av->tail;
    int v = av->head;
    int tUID = (*Tree::nodes)[u].memberOfTreeID;
    tree_type tUType = (*Tree::nodes)[u].memberOfTreeType;
    int tVID = (*Tree::nodes)[v].memberOfTreeID;
    tree_type tVType = (*Tree::nodes)[v].memberOfTreeType;

    // DEBUG:
//    printf("Arc being removed is: (%d,%d)\n", u, v);
//    fflush(stdout);

    if ((tUID != tVID) || (tUType != tVType)) {
        printf("Error: Attempting to remove arc spanning two trees\n");
        exit(-1);
    }

    Tree *tUV = NULL;
    if (tUType == TYPE_I) {
        tUV = (*setBtI)[tUID];
    } else { // tUType == TYPE_II
        tUV = (*setBtII)[tUID];
    }
  
    if ((tUV->treeType == TYPE_I) && 
        (av->varIndex == tUV->extraArc->varIndex)) {
        // Remove extra arc from type I tree and convert to type II [O(n)]
        Tree::convertTItoTII(tUV);
    } else { 
        // Things are a bit more complicated. Split the tree into two type II 
        // trees, then add the extra arc (if needed) to either create a type I 
        // tree and a type II tree or to merge the two type II trees into a 
        // type II tree. [O(|Ti| + |Tj| + n) = O(n)]
        Tree::removeTree(tUV);
        Tree::splitTree(tUV, av);
        if (tUV->treeType == TYPE_I) {
            Tree::insertIntoTrees(tUV->extraArc);
        }
        delete tUV;
    }

    return;
}

// Running time: O(|Ti|)
void Tree::splitTree(Tree *tX, ArcVar *removedArc)
{
    TreeTypeII *tTop = Tree::createNewTree(TYPE_II);
    TreeTypeII *tBottom = Tree::createNewTree(TYPE_II);
    TreeTypeII *tCur = tTop;

    stack<int> curNodeStack;
    stack<int> prevNodeStack;

    curNodeStack.push(tX->root);
    prevNodeStack.push(-1);

    while (!curNodeStack.empty()) {
        int curNode = curNodeStack.top(); curNodeStack.pop();
        if (curNode < 0) { // We've finished processing the bottom tree
            tCur = tTop;
            if (curNodeStack.empty()) {
                break;
            }
            curNode = curNodeStack.top(); curNodeStack.pop();
        }
        int prevNode = prevNodeStack.top(); prevNodeStack.pop();

        Node &tnC = (*Tree::nodes)[curNode];
        tnC.memberOfTreeID = tCur->treeID;
        tnC.memberOfTreeType = tCur->treeType;
        if (tCur->nodesInTree.size() == 0) {
            tCur->root = curNode; 
        }
        tCur->nodesInTree.push_back(curNode);

        int removedArcInd = -1;
        for (int arcInd = 0; arcInd < tnC.incidentArcs.size(); ++arcInd) {
            ArcVar *nextArc = tnC.incidentArcs[arcInd];
            if (nextArc->varIndex == removedArc->varIndex) {
                removedArcInd = arcInd;
                continue;
            }
            int nextNode = (nextArc->tail == curNode) ? nextArc->head 
                                                      : nextArc->tail;
            if (nextNode != prevNode) {
                tCur->arcsInTree.push_back(nextArc);
                curNodeStack.push(nextNode);
                prevNodeStack.push(curNode);
            }
        }

        // Now we check if we need to delete the removed arc from this node's 
        // incidence list
        if (removedArcInd >= 0) {
            tnC.incidentArcs[removedArcInd] = tnC.incidentArcs.back();
            tnC.incidentArcs.pop_back();
            // If we haven't started the bottom tree, do so now
            if (prevNode > -2) { // Switch to processing the bottom tree
                int nextNode = (removedArc->tail == curNode) 
                                        ? removedArc->head : removedArc->tail;
                curNodeStack.push(-1);
                curNodeStack.push(nextNode);
                prevNodeStack.push(-2);
                tCur = tBottom;
            }
        }
    }
    return;
}

/*****************************************************************************/
/* Update Trees Methods                                                      */
/*****************************************************************************/
// Running time: O(n)
void Tree::updateTrees()
{
    for (int tIind = 0; tIind < setBtI->size(); ++tIind) {
        TreeTypeI *tI = (*setBtI)[tIind];
        if (!tI->upToDate) { // [O(|Ti|)]
            // Set root of tree to tail of extra arc to keep root on cycle
            tI->root = tI->extraArc->tail;
            tI->initializeTreeIndices();
            tI->upToDate = true;
        }
    }
    for (int tIIind = 0; tIIind < setBtII->size(); ++tIIind) {
        TreeTypeI *tII = (*setBtII)[tIIind];
        if (!tII->upToDate) { // [O(|Ti|)]
            tII->initializeTreeIndices();
            tII->upToDate = true;
        }
    }
    return;
}

/*****************************************************************************/
/* Dynamic tree methods needed for computing potentials and flows            */
/*****************************************************************************/
// Running time: O(|Ti|)
void Tree::initializeTreeIndices()
{
    // Initialize depth, pred, and thread indices
    stack<int> curNodeStack;
    stack<ArcVar *> predArcStack;

    threads.clear(); 
    curNodeStack.push(root);

    while (!curNodeStack.empty()) {
        int curNode = curNodeStack.top(); curNodeStack.pop();
        Node &tnC = (*Tree::nodes)[curNode];

        int predNode = -1;
        if (threads.size() > 0) {
            ArcVar *predArc = predArcStack.top(); predArcStack.pop();
            predNode = (curNode == predArc->tail) ? predArc->head 
                                                  : predArc->tail;
            tnC.depth = (*Tree::nodes)[predNode].depth + 1;
            tnC.predArc = predArc;
            (*Tree::nodes)[threads.back()].thread = curNode;
        } else { // At root
            tnC.depth = 0;
            tnC.predArc = NULL;     
        }
        tnC.pred = predNode;
        threads.push_back(curNode);

        for (int arcInd = tnC.incidentArcs.size()-1; arcInd >= 0; --arcInd) {
            ArcVar *av = tnC.incidentArcs[arcInd];
            int u = av->tail;
            int v = av->head;
            if ((u == predNode) || (v == predNode)) {
                continue;
            }
            if (curNode == u) {
                curNodeStack.push(v);
            } else { // curNode == v
                curNodeStack.push(u);
            }
            predArcStack.push(av);
        }
    }

    (*Tree::nodes)[threads.back()].thread = root;

    return;
}

