// Copyright (c) 2022 Yi Pan <ypan1988@gmail.com>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  A copy of the GNU General Public License is available at
//  https://www.R-project.org/Licenses/

#ifndef FUNCTOR_H_
#define FUNCTOR_H_

#include <RcppArmadillo.h>

namespace rminqa {
class Functor {
public:
  Functor() {}
  virtual ~Functor() {}
  virtual double operator()(const arma::vec &par) = 0;

  int feval = 0;
};

inline double minqa_objfun(long n, const double *x, void *data) {
  arma::vec par(x, n);
  ++(static_cast<Functor *>(data)->feval);
  return static_cast<Functor *>(data)->operator()(par);
}

}

#endif
