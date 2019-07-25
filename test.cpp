#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <ctime>
//#include <omp.h>
#include <boost/math/special_functions/digamma.hpp>
#include "KdTree.h"

using namespace std;

int main(){
  /* prepare data here */
  double sample_N = 400000;
  int K = 100;
  vector<double> test(sample_N);
  vector<double> test2x(sample_N);
  vector<double> test2y(sample_N);
  vector<double> test3x(sample_N);
  vector<double> test3y(sample_N);
  vector<double> test3z(sample_N);

  random_device rnd;
  mt19937 mt(rnd());
  normal_distribution<> norm(0, 0.02);
  normal_distribution<> norm2(0, 1);

  double temp;
  for(int j=0; j<sample_N; j++){
    test[j] = norm(mt);
    
    temp = norm2(mt);
    test2x[j] = 2*temp;
    test2y[j] = temp + 2*norm2(mt);

    test3x[j] = norm2(mt);
    test3y[j] = norm2(mt);
    test3z[j] = norm2(mt);
  }

  string filename = "test.txt";
  ofstream ofs(filename);
  for(int j=0; j<sample_N; j++){
    ofs << setprecision(15) << test2x[j] << ' ' << test2y[j] << ' ' << test[j] << '\n';
  }
  ofs.close();


  /* construct Kd-tree here */
  clock_t start, end;
  start = clock();
  KdTree kdTree = KdTree(vector<vector<double> >{test});
  KdTree kdTree2 = KdTree(vector<vector<double> >{test2x, test2y});
  KdTree kdTree2x = KdTree(vector<vector<double> >{test2x});
  KdTree kdTree2y = KdTree(vector<vector<double> >{test2y});
  KdTree kdTree3 = KdTree(vector<vector<double> >{test3x, test3y, test3z});
  end = clock();
  cout << "construct" << endl;
  cout << "＊time = " << (double)(end-start)/CLOCKS_PER_SEC << "sec\n";

  
  /* 
     Basic usage of our library
     We use test3 (3d point cloud) here
  */
  double ret = 0;
  double eps = 0;
  int Nx = 0;
  int Ny = 0;

  //the distance to the Kth nearest neighbour from the first (index = 0) data point
  start = clock();
  {
    eps = kdTree3.findKNearestNeighbour(K, 0);
  }
  end = clock();
  cout << "find 20thNN" << endl;
  cout << "＊time = " << (double)(end-start)/CLOCKS_PER_SEC << "sec\n";

  //the number of data points which exist within r < 1.0 from the first(0) data point
  start = clock();
  {
    Nx = kdTree3.countPointsWithinR(0, 1);
  }
  end = clock();
  cout << "count points N = " << Nx << endl;
  cout << "＊time = " << (double)(end-start)/CLOCKS_PER_SEC << "sec\n";
  
  /*
    KL method 
    We calculate the Shannon entropy of test
  */
  start = clock();
  {
    #pragma omp parallel for reduction(+:ret) private(eps)
    for(int j=0; j<sample_N; j++){
      eps = kdTree.findKNearestNeighbour(K, j);
      ret += log(2*eps)/sample_N;
    }
    ret = ret - boost::math::digamma(K) + boost::math::digamma(sample_N);
  }
  end = clock();
  cout << "Shannon entropy " << ret << endl;
  cout << "＊analytical value = " << -2.493084472 << endl; //1/2 + log(sqrt(2*pi)*0.02)
  cout << "＊time = " << (double)(end-start)/CLOCKS_PER_SEC << "sec /thread_num\n";

  
  /*
    KSG method 
    We calculate the Mutual information between test2x and test2y
  */
  start = clock();
  ret = 0;
  {
    #pragma omp parallel for reduction(+:ret) private(eps, Nx, Ny)
    for(int j=0; j<sample_N; j++){
      eps = kdTree2.findKNearestNeighbour(K, j);
      Nx = kdTree2x.countPointsWithinR(j, eps);
      Ny = kdTree2y.countPointsWithinR(j, eps);
      ret -= boost::math::digamma(Nx+1)/sample_N + boost::math::digamma(Ny+1)/sample_N;
    }
    ret = ret + boost::math::digamma(K) + boost::math::digamma(sample_N);
  }
  end = clock();
  cout << "Mutual information " << ret << endl;
  cout << "＊analytical value = " << 0.111571775 << endl; //-1/2*log(4/5)
  cout << "＊time = " << (double)(end-start)/CLOCKS_PER_SEC << "sec /thread_num\n";
  
  return 0;
}
