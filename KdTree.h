#ifndef _INC_KDTREE
#define _INC_KDTREE
#include "KdTree.h"
using namespace std;

namespace NormType
{
  static const int MAX = 0;
  static const int EUCLIDIAN = 1;
}


class KdTreeNode
{
 public:
  int index;
  KdTreeNode* leftTree;
  KdTreeNode* rightTree;
  KdTreeNode() {}
  KdTreeNode(int, KdTreeNode*, KdTreeNode*);
};


class KdTree
{
 private:
  int totalDimensions;
  int numObservations;
  int rootNodeIndex;
  int normType;
  int K;
  KdTreeNode endOfTree;
  vector<vector<double> > data;
  vector<vector<double> > dataForNormCalc;
  vector<KdTreeNode> tree;
 
 public:
  KdTree(vector<vector<double> >);
  int constructKdTree(int, int, vector<vector<int> >&);
  void setNormType(int n) { normType = n; }
  double norm(int, int);
  double findKNearestNeighbour(int, int);
  void findKNearestNeighboursCalc(int, vector<pair<double, int> >&, double&, int&, KdTreeNode&, int);
  void fastSort(vector<pair<double, int> >&);
  int countPointsWithinR(int, double);
  void countPointsWithinRCalc(int, double, KdTreeNode&, int&, int);
  //void countPointsWithinRUsingCResult();
};
#endif
