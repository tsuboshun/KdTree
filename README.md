# KdTree for calculating Shannon entropy
## 1. About
This is a simple code for the calculation of Shannon entropy or mutual information, written in C++.
You can implement [KL method](http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=ppi&paperid=797&option_lang=eng) for the calculation of Shannon entropy and [KSG method](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.69.066138) for that of mutual information.
Please refer to test.cpp for the usage.
I note that there is already a famous library called [JIDT](https://github.com/jlizier/jidt) for these methods, which is written in Java.
I wrote this code referring to the code of JIDT, and checked the consistency of the calculation results.
This code can compute Shannon entropy faster than JIDT, but it is around six times slower than JIDT for the calculation of mutual information.

## 2. How to run
You need to install the [Boost C++ library](https://www.boost.org) beforehand. If you are in the directory of KdTree.cpp, you can compile test.cpp as follows.  
```g++ test.cpp KdTree.cpp -I/path_to_boost_library -std=c++11```
