#ifndef _INC_KDTREE
#define _INC_KDTREE
#include "KdTree.h"
using namespace std;

struct KdTreeNode
{
  int index;
  KdTreeNode* leftTree;
  KdTreeNode* rightTree;
};

class KdTree
{
 private:
  int totalDimensions;
  int numObservations;
  int rootNodeIndex;
  int K;
  KdTreeNode endOfTree;
  vector<vector<double> > data;
  vector<vector<double> > dataForNormCalc;
  struct KdTreeNode* tree;
 
 public:
  KdTree(vector<vector<double> >);
  ~KdTree();
  int constructKdTree(int, int&, vector<vector<int> >&);
  double norm(int&, int&);
  double findKNearestNeighbour(int, int);
  void findKNearestNeighboursCalc(int&, vector<pair<double, int> >&, double&, int&, KdTreeNode&, int);
  void fastSort(vector<pair<double, int> >&);
  int countPointsWithinR(int, double);
  void countPointsWithinRCalc(int&, double&, KdTreeNode&, int&, int);
};
#endif
