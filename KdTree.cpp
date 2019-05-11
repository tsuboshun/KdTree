#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "KdTree.h"
using namespace std;


inline bool compare_pair(pair<double, int>& b1, pair<double, int>& b2){
  return b1.first > b2.first || ( b1.first == b2.first && b1.second > b2.second);
}

KdTree::KdTree(vector<vector<double> > d)
{
  data = d;
  totalDimensions = data.size();
  numObservations = data[0].size();
  dataForNormCalc = vector<vector<double> >(numObservations, vector<double>(totalDimensions, 0));
  for(int j=0; j<totalDimensions; ++j){
    for(int i=0; i<numObservations; ++i){
      dataForNormCalc[i][j] = data[j][i];
    }
  }
  vector<vector<int> > masterSortedArrayIndices 
    = vector<vector<int> >(totalDimensions, vector<int>(numObservations, 0));
  vector<pair<double, int> > distTable = 
    vector<pair<double, int> >(numObservations, make_pair(0., 0));
  for(int i=0; i<totalDimensions; ++i){
    for(int j=0; j<numObservations; ++j){
      distTable[j].second = j;
      distTable[j].first = data[i][j];
    }
    sort(distTable.begin(), distTable.end());
    for(int j=0; j<numObservations; ++j){
      masterSortedArrayIndices[i][j] = distTable[j].second;
    }
  }
  tree = (struct KdTreeNode*) malloc(sizeof(KdTreeNode)*numObservations);
  for(int j=0; j<numObservations; ++j)
    tree[j].index = j;
  endOfTree = {-1, NULL, NULL};
  rootNodeIndex = constructKdTree(0, numObservations, masterSortedArrayIndices);
}


KdTree::~KdTree()
{
  vector<vector<double> >().swap(data);
  vector<vector<double> >().swap(dataForNormCalc);
  free(tree);
}


inline int KdTree::constructKdTree(int currentDim, int& numPoints, vector<vector<int> >& sortedArrayIndices)
{
  int temp;
  if(numPoints == 1){
    temp = sortedArrayIndices[currentDim][0];
    tree[temp].leftTree = &endOfTree;
    tree[temp].rightTree = &endOfTree;
    return temp;
  }
  
  if(numPoints == 2){
    temp = sortedArrayIndices[currentDim][1];
    tree[temp].leftTree = &endOfTree;
    tree[temp].rightTree = &endOfTree;
    temp = sortedArrayIndices[currentDim][0];
    tree[temp].leftTree = &endOfTree;
    tree[temp].rightTree = &tree[sortedArrayIndices[currentDim][1]];
    return temp;
  }
  
  int candidateSplitPoint = numPoints/2;
  while((candidateSplitPoint > 0) &&
	(data[currentDim][sortedArrayIndices[currentDim][candidateSplitPoint-1]] ==
	 data[currentDim][sortedArrayIndices[currentDim][candidateSplitPoint]])){
    --candidateSplitPoint;
  }
  double medianValueInCurrentDim = data[currentDim][sortedArrayIndices[currentDim][candidateSplitPoint]];
  int sampleNumberForSplitPoint = sortedArrayIndices[currentDim][candidateSplitPoint];
  int leftNumPoints = candidateSplitPoint;
  int rightNumPoints = numPoints-candidateSplitPoint-1;

  vector<vector<int> > newSortedArrayIndicesLeft = vector<vector<int> >
    (totalDimensions, vector<int>(leftNumPoints, 0));
  vector<vector<int> > newSortedArrayIndicesRight = vector<vector<int> >
    (totalDimensions, vector<int>(rightNumPoints, 0));

  int leftIndex;
  int rightIndex;
  int sampleNumberInData;
  for(int dim=0; dim<totalDimensions; dim++){
    if(dim == currentDim){
      for(int i=0; i<leftNumPoints; ++i)
	newSortedArrayIndicesLeft[dim][i] = sortedArrayIndices[dim][i];
      for(int i=0; i<rightNumPoints; ++i)
	newSortedArrayIndicesRight[dim][i] = sortedArrayIndices[dim][candidateSplitPoint+1+i];
      continue;
    }
    leftIndex = 0;
    rightIndex = 0;
    for(int i=0; i<numPoints; ++i){
      sampleNumberInData = sortedArrayIndices[dim][i];
      if (sampleNumberInData == sampleNumberForSplitPoint){
	continue;
      }
      if(data[currentDim][sampleNumberInData] < medianValueInCurrentDim){
	newSortedArrayIndicesLeft[dim][leftIndex++] = sampleNumberInData;
      } else {
	newSortedArrayIndicesRight[dim][rightIndex++] = sampleNumberInData;
      }
    }
  }
  int newDim = (currentDim + 1) % totalDimensions;
  
  if(leftNumPoints==0)
    tree[sampleNumberForSplitPoint].leftTree = &endOfTree;
  else
    tree[sampleNumberForSplitPoint].leftTree = &tree[constructKdTree(newDim, leftNumPoints, newSortedArrayIndicesLeft)];
  tree[sampleNumberForSplitPoint].rightTree = &tree[constructKdTree(newDim, rightNumPoints, newSortedArrayIndicesRight)];
  return sampleNumberForSplitPoint;
}


inline double KdTree::norm(int& i1, int& i2)
{
  double dist = 0.;
  double diff = 0.;
  for(int d=0; d<totalDimensions; d++){
    diff = abs(dataForNormCalc[i1][d] - dataForNormCalc[i2][d]);
    if(diff > dist)
      dist = diff;
  }
  return dist;
}


double KdTree::findKNearestNeighbour(int k, int sampleIndex)
{
  K = k;
  vector<pair<double, int> > currentKBestVec = vector<pair<double, int> >(K, pair<double, int>(-1, -1));
  double currentKBestDist = 0.;
  int currentKBestIndex = 0;
  findKNearestNeighboursCalc(sampleIndex, currentKBestVec, currentKBestDist, currentKBestIndex, tree[rootNodeIndex], 0);
  return currentKBestVec[0].first;
}


//currentKBest is to be sorted in a descending order
inline void KdTree::findKNearestNeighboursCalc(int& sampleIndex, vector<pair<double, int> >& currentKBestVec, double& currentKBestDist, int& currentKBestIndex, KdTreeNode& node, int level)
{
  if(node.index == -1)
    return;
  
  double Dist = norm(sampleIndex, node.index);
  if(node.index!=sampleIndex){ //pass if it is the target data itself
    if(currentKBestIndex < K){
      currentKBestVec[currentKBestIndex].first = Dist;
      currentKBestVec[currentKBestIndex].second = node.index;
      currentKBestIndex++;
      if(currentKBestIndex==K){
	sort(currentKBestVec.begin(), currentKBestVec.end(), compare_pair);
	currentKBestDist = currentKBestVec[0].first;
      }
    }else if(currentKBestDist > Dist){
      currentKBestVec[0].first = Dist;
      currentKBestVec[0].second = node.index;
      fastSort(currentKBestVec);
      currentKBestDist = currentKBestVec[0].first;
    }
  }
  
  int nextLevel = (level + 1)%totalDimensions;
  double dist = data[level][sampleIndex]-data[level][node.index];
  if(dist<=0){
    findKNearestNeighboursCalc(sampleIndex, currentKBestVec, currentKBestDist, currentKBestIndex, *node.leftTree, nextLevel);
    if(currentKBestIndex < K || currentKBestDist > -dist)
      findKNearestNeighboursCalc(sampleIndex, currentKBestVec, currentKBestDist, currentKBestIndex, *node.rightTree, nextLevel);
  }else{
    findKNearestNeighboursCalc(sampleIndex, currentKBestVec, currentKBestDist, currentKBestIndex, *node.rightTree, nextLevel);
    if(currentKBestIndex < K || dist < currentKBestDist)
      findKNearestNeighboursCalc(sampleIndex, currentKBestVec, currentKBestDist, currentKBestIndex, *node.leftTree, nextLevel);
  }
  return;
}


inline void KdTree::fastSort(vector<pair<double, int> >& vec)
{ 
  pair<double, int> vec0 = vec[0];
  for(int i=0; i<vec.size()-1; ++i){
    if(vec0.first >= vec[i+1].first){
      vec[i] = vec0;
      return;
    }
    vec[i] = vec[i+1];
  }
  vec.back() = vec0;
  return;
}

  
int KdTree::countPointsWithinR(int sampleIndex, double R) //how many points there are within r < R from the target data (sampleIndex) excluding itself
{
  int count = 0;
  countPointsWithinRCalc(sampleIndex, R, tree[rootNodeIndex], count, 0);
  return count-1; //excluding the target data
}


inline void KdTree::countPointsWithinRCalc(int& sampleIndex, double& R, KdTreeNode& node, int& count, int level)
{
  if(node.index==-1)
    return;
  
  for(int d=0; d<totalDimensions; d++){
    if(abs(dataForNormCalc[sampleIndex][d]-dataForNormCalc[node.index][d]) >= R)
      break;
    if(d==totalDimensions-1)
      count++;
  }
  
  double dist = dataForNormCalc[sampleIndex][level]-dataForNormCalc[node.index][level];
  int nextLevel = (level+1) % totalDimensions;
  
  if(dist<=0){
    countPointsWithinRCalc(sampleIndex, R, *node.leftTree, count, nextLevel);
    if(-dist < R)
      countPointsWithinRCalc(sampleIndex, R, *node.rightTree, count, nextLevel);
  }else{
    countPointsWithinRCalc(sampleIndex, R, *node.rightTree, count, nextLevel);
    if(dist < R)
      countPointsWithinRCalc(sampleIndex, R, *node.leftTree, count, nextLevel);
  }  
  return;
}
