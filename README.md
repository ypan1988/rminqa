rminqa: Derivative-Free Optimization in R using C++.
====

[![Build Status](https://travis-ci.org/ypan1988/rminqa.svg?branch=master)](https://travis-ci.org/ypan1988/rminqa)
[![cran version](http://www.r-pkg.org/badges/version/rminqa)](https://cran.r-project.org/web/packages/rminqa)
[![downloads](http://cranlogs.r-pkg.org/badges/rminqa)](http://cranlogs.r-pkg.org/badges/rminqa)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/rminqa)](http://cranlogs.r-pkg.org/badges/grand-total/rminqa)

## Features

* Perform derivative-free optimization algorithms in R using C++.
* A wrapper interface is provided to call `C` function of the 'bobyqa' implementation (See <https://github.com/emmt/Algorithms/tree/master/bobyqa>).

## Installation

Get the development version from github:
```R
install.packages("devtools")
library(devtools)
devtools::install_github("ypan1988/rminqa", dependencies=TRUE)
```

Or the released version from CRAN:
```R
install.packages("rminqa")
```

## A Quick Example
An example of using `minqa::bobyqa()` in `R` environment:
```R
fr <- function(x) {   ## Rosenbrock Banana function
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
(x1 <- minqa::bobyqa(c(1, 2), fr, lower = c(0, 0), upper = c(4, 4)))
## => optimum at c(1, 1) with fval = 0
str(x1) # see that the error code and msg are returned
```

Corresponding code written in `C++` using package `rminqa` (file `demo.cpp`):
```cpp
#include <cmath>  // std::pow

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "rbobyqa.h"
// [[Rcpp::depends(rminqa)]]

using namespace rminqa;

class Rosen : public Functor {
public:
  double operator()(const arma::vec &x) override {
    double x1 = x(0);
    double x2 = x(1);
    return 100 * std::pow((x2 - x1 * x1), 2) + std::pow(1 - x1, 2);
  }
};

// [[Rcpp::export]]
void bobyqa_rosen() {
  Rosen rb;
  Rbobyqa<Rosen> opt;

  arma::vec x = {1, 2};
  opt.minimize(rb, x);

  Rcpp::Rcout << "-------------------------" << std::endl;
  Rcpp::Rcout << "par = \n" << opt.par() << std::endl;
  Rcpp::Rcout << "fval = " << opt.fval() << std::endl;
  Rcpp::Rcout << "feval = " << opt.feval() << std::endl;
  Rcpp::Rcout << "msg = " << opt.msg() << std::endl;
  Rcpp::Rcout << "-------------------------" << std::endl;
}
```

Compile and run the function in `R`:
```R
library(Rcpp)
sourceCpp("~/demo.cpp") # you may need to change the directory
bobyqa_rosen()
```

Then you will get expected output as follows:
```
-------------------------
par = 
   1.0000
   1.0000

fval = 9.98796e-15
feval = 341
msg = Normal exit from bobyqa
-------------------------
```
