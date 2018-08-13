#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "KdTree.h"
using namespace std;

bool compare_pair(pair<double, int> b1, pair<double, int> b2){
  if(b1.first > b2.first)
    return true;
  else if(b1.first < b2.first)
    return false;
  else
    if(b1.second > b2.second)
      return true;
    else
      return false;
}

KdTreeNode::KdTreeNode(int i, KdTreeNode* lt, KdTreeNode* rt)
{
    index = i;
    leftTree = lt;
    rightTree = rt;
}

KdTree::KdTree(vector<vector<double> > d)
{
  normType = NormType::MAX;
  data = d;
  totalDimensions = data.size();
  numObservations = data[0].size();
  dataForNormCalc = vector<vector<double> >(numObservations, vector<double>(totalDimensions, 0));
  for(int i=0; i<numObservations; i++){
    for(int j=0; j<totalDimensions; j++){
      dataForNormCalc[i][j] = data[j][i];
    }
  }
  vector<vector<int> > masterSortedArrayIndices 
    = vector<vector<int> >(totalDimensions, vector<int>(numObservations, 0));
  vector<pair<double, int> > distTable = 
    vector<pair<double, int> >(numObservations, make_pair(0., 0));
  for(int i=0; i<totalDimensions; i++){
    for(int j=0; j<numObservations; j++){
      distTable[j].second = j;
      distTable[j].first = data[i][j];
    }
    sort(distTable.begin(), distTable.end());
    for(int j=0; j<numObservations; j++){
      masterSortedArrayIndices[i][j] = distTable[j].second;
    }
  }
  tree.resize(numObservations);
  for(int j=0; j<numObservations; j++)
    tree[j].index = j;
  endOfTree = KdTreeNode(-1, NULL, NULL);
  rootNodeIndex = constructKdTree(0, numObservations, masterSortedArrayIndices);
}


int KdTree::constructKdTree(int currentDim, int numPoints, vector<vector<int> >& sortedArrayIndices)
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
    candidateSplitPoint--;
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
      for(int i=0; i<leftNumPoints; i++)
	newSortedArrayIndicesLeft[dim][i] = sortedArrayIndices[dim][i];
      for(int i=0; i<rightNumPoints; i++)
	newSortedArrayIndicesRight[dim][i] = sortedArrayIndices[dim][candidateSplitPoint+1+i];
      continue;
    }
    leftIndex = 0;
    rightIndex = 0;
    for(int i=0; i<numPoints; i++){
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

 double KdTree::norm(int i1, int i2)
{
  vector<double> x1 = dataForNormCalc[i1];
  vector<double> x2 = dataForNormCalc[i2];
  double dist = 0.;
  double diff = 0.;
  switch(normType){
  case NormType::MAX:
    for(int d=0; d<totalDimensions; d++){
      diff = x1[d] - x2[d];
      if(diff<0)
	diff = -diff;
      if(diff > dist)
	dist = diff;
    }
    return dist;
  case NormType::EUCLIDIAN: //TODO EUCLIDEANだとKd木の探索領域をもう少し減らせる
    for(int d=0; d<x1.size(); d++){
      diff = x1[d]-x2[d];
      dist += diff * diff;
    }
    return dist;
  }
}

//並列化できるように書くこと
double KdTree::findKNearestNeighbour(int k, int sampleIndex)
{
  K = k;
  vector<pair<double, int> > currentKBestVec = vector<pair<double, int> >(K, pair<double, int>(-1, -1));
  double currentKBestDist = 0.;
  int currentKBestIndex = 0;
  findKNearestNeighboursCalc(sampleIndex, currentKBestVec, currentKBestDist, currentKBestIndex, tree[rootNodeIndex], 0);
  return currentKBestVec[0].first;
}

//currentKBestは降順にソートする
void KdTree::findKNearestNeighboursCalc(int sampleIndex, vector<pair<double, int> >& currentKBestVec, double& currentKBestDist, int& currentKBestIndex, KdTreeNode& node, int level)
{
  if(node.index == -1)
    return;
  
  double Dist = norm(sampleIndex, node.index);
  if(Dist!=0){ //自分自身の場合はpass
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
  }else{
    findKNearestNeighboursCalc(sampleIndex, currentKBestVec, currentKBestDist, currentKBestIndex, *node.rightTree, nextLevel);
  }

  if(normType == NormType::EUCLIDIAN){
    if(currentKBestIndex < K || dist*dist < currentKBestDist){
      if(dist <= 0)
	findKNearestNeighboursCalc(sampleIndex, currentKBestVec, currentKBestDist, currentKBestIndex, *node.rightTree, nextLevel);
      else
	findKNearestNeighboursCalc(sampleIndex, currentKBestVec, currentKBestDist, currentKBestIndex, *node.leftTree, nextLevel);
    }
  }else{
    if(dist <= 0 && (currentKBestIndex < K || currentKBestDist > -dist))
      findKNearestNeighboursCalc(sampleIndex, currentKBestVec, currentKBestDist, currentKBestIndex, *node.rightTree, nextLevel);
    else if(dist > 0 && (currentKBestIndex < K || dist < currentKBestDist))
      findKNearestNeighboursCalc(sampleIndex, currentKBestVec, currentKBestDist, currentKBestIndex, *node.leftTree, nextLevel);
  }
  return;
}

void KdTree::fastSort(vector<pair<double, int> >& vec)
{ 
  pair<double, int> vec0 = vec[0];
  if(vec0.first > vec[1].first)
    return;
  if(vec0.first < vec.back().first){
    for(int i=0; i<vec.size()-1; i++)
      vec[i] = vec[i+1];
    vec.back() = vec0;
    return;
  }

  for(int i=2; i<vec.size(); i++){
    if(vec0.first >= vec[i].first){
      for(int j=0; j<i-1; j++)
	vec[j] = vec[j+1];
      vec[i-1]=vec0;
      return;
    }
  }
}


int KdTree::countPointsWithinR(int sampleIndex, double R) //自分自信を除いてR以内に何点あるか
{
  int count = 0;
  countPointsWithinRCalc(sampleIndex, R, tree[rootNodeIndex], count, 0);
  return count-1; //最後に自分自身を引く
}

void KdTree::countPointsWithinRCalc(int sampleIndex, double R, KdTreeNode& node, int& count, int level)
{
  if(node.index==-1)
    return;
  
  vector<double> x1 = dataForNormCalc[sampleIndex];
  vector<double> x2 = dataForNormCalc[node.index];
  double Dist;
  for(int d=0; d<totalDimensions; d++){
    Dist = x1[d]-x2[d];
    if(Dist < 0)
      Dist = -Dist;
    if(Dist >= R)
      break;
    else if(d==totalDimensions-1)
      count++;
  }//EUCLIDEANの場合を後で作る

  int nextLevel = (level+1) % totalDimensions;
  double dist = x1[level]-x2[level];
  
  if(dist<=0)
    countPointsWithinRCalc(sampleIndex, R, *node.leftTree, count, nextLevel);
  else
    countPointsWithinRCalc(sampleIndex, R, *node.rightTree, count, nextLevel);
  
  if(normType == NormType::EUCLIDIAN){
    if(dist*dist < R){
      if(dist <= 0)
	countPointsWithinRCalc(sampleIndex, R, *node.rightTree, count, nextLevel);
      else
	countPointsWithinRCalc(sampleIndex, R, *node.leftTree, count, nextLevel);
    }
  }else{
    if(dist <= 0 && -dist < R)
      countPointsWithinRCalc(sampleIndex, R, *node.rightTree, count, nextLevel);
    else if(dist > 0 && dist < R)
      countPointsWithinRCalc(sampleIndex, R, *node.leftTree, count, nextLevel);
  }
  return;
}
