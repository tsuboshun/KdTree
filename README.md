# KdTree for calculating Shannon entropy
## 1. About
It is simple and an easy to use library for calculating Shannon entropy or mutual information, written with C++.
You can implement [KL method](http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=ppi&paperid=797&option_lang=eng) for the calculation of Shannon entropy and [KSG method](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.69.066138) for that of mutual information.
Please refer to test.cpp for the usage.
As similar libraries, there are [nanoflann](https://github.com/jlblancoc/nanoflann) which is a C++ header-only library but may require you to write an adapter class or [JIDT](https://github.com/jlizier/jidt) which is a Java library.
I wrote this library referring to the code of JIDT, and checked the consistency of the calculation results.
My code is faster than JIDT in terms of the calculation of Shannon entropy, but around three times slower in terms of the KSG method right now. I welcome any pull requests.

## 2. How to run
You need to install the Boost C++ library beforehand. If you are in the directory of KdTree.cpp, you can compile test.cpp as follows.  
```g++ test.cpp KdTree.cpp -I/path_to_boost_library ```
