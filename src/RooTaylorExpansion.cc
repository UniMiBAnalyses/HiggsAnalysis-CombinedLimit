#include "RooFit.h"

#include "Riostream.h"
#include <math.h>

#include "../interface/RooTaylorExpansion.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooMath.h"
#include "RooConstVar.h"
#include "RooListProxy.h"

using namespace std;

ClassImp(RooTaylorExpansion)
  ;

//_____________________________________________________________________________
RooTaylorExpansion::RooTaylorExpansion(const char* name, const char* title,
                                       const RooArgList& x,
                                       const RooArgList& x0, unsigned int order,
                                       const TVectorD& terms)
    : RooAbsReal(name, title),
      _x("x", "Observables", this, kTRUE, kFALSE),
      _x0("x0", "initial point", this, kTRUE, kFALSE),
      _order(order),
      _terms(terms) {
  _x.add(x);
  _x0.add(x0);
}

//_____________________________________________________________________________
RooTaylorExpansion::RooTaylorExpansion(const RooTaylorExpansion& other,
                                       const char* name)
    : RooAbsReal(other, name),
      _x("x", this, other._x),
      _x0("x0", this, other._x0),
      _order(other._order),
      _terms(other._terms) {}

//_____________________________________________________________________________
Double_t RooTaylorExpansion::evaluate() const {
  double fullsum = 0.;
  unsigned pos = 0;
  unsigned int factorial = 1;

  unsigned nx = _x.getSize();
  std::vector<double> dx(nx);

  for (unsigned i = 0; i < nx; ++i) {
    dx[i] = ((RooAbsReal*)(_x.at(i)))->getVal() - ((RooAbsReal*)(_x0.at(i)))->getVal();
  }
  for (unsigned n = 0; n <= _order; ++n) {
    double order_term = 1. / double(factorial);
    factorial *= (n+1);

    std::vector<unsigned> tracker(n, 0);
    double order_sum = 0.;
    for (unsigned i = 0; i < unsigned(std::pow(nx, n) + 0.5); ++i) {
      double term = _terms[pos];
      for (unsigned t = 0; t < n; ++t) {
        term *= dx[tracker[t]];
      }
      order_sum += term;

      // increment the rightmost digit of the tracker
      for (unsigned t = 0; t < n; ++t) {
        if (tracker[n - 1 - t] == (nx - 1)) {
          tracker[n - 1 - t] = 0;
        } else {
          ++tracker[n - 1 - t];
          break;
        }
      }
      ++pos;
    }
    fullsum += (order_term * order_sum);
  }
  return fullsum;
}
