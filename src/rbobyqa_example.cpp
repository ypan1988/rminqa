#include "RcppArmadillo.h"

#include "rbobyqa.h"
using namespace rminqa;

class Rosen : public Functor {
public:
  double operator()(const arma::vec &x) override {
    double x1 = x(0);
    double x2 = x(1);
    return 100 * std::pow((x2 - x1 * x1), 2) + std::pow(1 - x1, 2);
  }
};

//'@title Example 1a: Minimize Rosenbrock function using bobyqa
//'@description Minimize Rosenbrock function using bobyqa and
//'             expect a normal exit from bobyqa.
//'@return No return value, called for side effects.
//'@examples
//'fr <- function(x) {   ## Rosenbrock Banana function
//'  x1 <- x[1]
//'  x2 <- x[2]
//'  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
//'}
//'(x1 <- minqa::bobyqa(c(1, 2), fr, lower = c(0, 0), upper = c(4, 4)))
//'## => optimum at c(1, 1) with fval = 0
//'str(x1) # see that the error code and msg are returned
//'
//'## corresponding C++ implementation:
//'bobyqa_rosen_x1()
//'@export
// [[Rcpp::export]]
void bobyqa_rosen_x1() {
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

//'@title Example 1b: Minimize Rosenbrock function using bobyqa
//'@description Minimize Rosenbrock function using bobyqa and
//'             expect a normal exit from bobyqa.
//'@return No return value, called for side effects.
//'@examples
//'fr <- function(x) {   ## Rosenbrock Banana function
//'  x1 <- x[1]
//'  x2 <- x[2]
//'  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
//'}
//'# check the error exits
//'# too many iterations
//'x1e <- minqa::bobyqa(c(1, 2), fr, lower = c(0, 0), upper = c(4, 4), control = list(maxfun=50))
//'str(x1e)
//'
//'## corresponding C++ implementation:
//'bobyqa_rosen_x1e()
//'@export
// [[Rcpp::export]]
void bobyqa_rosen_x1e() {
  Rosen rb;
  Rbobyqa<Rosen> opt;
  opt.set_lower({0, 0});
  opt.set_upper({4, 4});
  opt.control.maxfun = 50;

  arma::vec x = {1, 2};
  opt.minimize(rb, x);

  Rcpp::Rcout << "-------------------------" << std::endl;
  Rcpp::Rcout << "par = \n" << opt.par() << std::endl;
  Rcpp::Rcout << "fval = " << opt.fval() << std::endl;
  Rcpp::Rcout << "feval = " << opt.feval() << std::endl;
  Rcpp::Rcout << "msg = " << opt.msg() << std::endl;
  Rcpp::Rcout << "-------------------------" << std::endl;
}