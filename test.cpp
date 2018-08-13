#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <boost/math/special_functions/digamma.hpp>
#include "KdTree.h"

using namespace std;

int main(){
  double sample_N = 100000;
  random_device rnd;
  mt19937 mt(rnd());
  normal_distribution<> norm(0, 0.02);
  normal_distribution<> norm2(0, 1);
  vector<double> test(sample_N);
  vector<double> test2x(sample_N);
  vector<double> test2y(sample_N);
  double temp;
  for(int j=0; j<sample_N; j++){
    test[j] = norm(mt);
    temp = norm2(mt);
    test2x[j] = 2*temp;
    test2y[j] = temp + 2*norm2(mt);
  }

  string filename = "test.txt";
  ofstream ofs(filename);

  for(int j=0; j<sample_N; j++){
    ofs << test2x[j] << ' ' << test2y[j] << '\n';
  }
  ofs.close();

  int K = 20;
  KdTree kdTree = KdTree(vector<vector<double> >{test});
  double ret = 0;
  double eps = 0;
  
  KdTree kdTree2 = KdTree(vector<vector<double> >{test2x, test2y});
  KdTree kdTree2x = KdTree(vector<vector<double> >{test2x});
  KdTree kdTree2y = KdTree(vector<vector<double> >{test2y});
  double ret2 = 0;
  double eps2 = 0;
  int Nx = 0;
  int Ny = 0;
  
  for(int j=0; j<sample_N; j++){
    eps = kdTree.findKNearestNeighbour(K, j);
    ret += log(2*eps)/sample_N;

    eps2 = kdTree2.findKNearestNeighbour(K, j);
    Nx = kdTree2x.countPointsWithinR(j, eps2);
    Ny = kdTree2y.countPointsWithinR(j, eps2);
    ret2 -= boost::math::digamma(Nx+1)/sample_N + boost::math::digamma(Ny+1)/sample_N;
  }
  ret = ret - boost::math::digamma(K) + boost::math::digamma(sample_N);
  ret2 = ret2 + boost::math::digamma(K) + boost::math::digamma(sample_N);
  cout << ret << " expected value: " << -2.493084472 << endl; //理論値 1/2 + log(sqrt(2*pi)*0.02) = -2.493084472
  cout << ret2 << " expected value: " << 0.111571775 << endl; //理論値 -1/2*log(4/5) = 0.111571775
  
  return 0;
}
