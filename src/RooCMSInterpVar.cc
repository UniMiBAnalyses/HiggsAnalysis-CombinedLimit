#include "HiggsAnalysis/CombinedLimit/interface/RooCMSInterpVar.h"

#include <fstream>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include "TMath.h"

RooCMSInterpVar::RooCMSInterpVar(const char *name, const char *title, RooAbsReal &xvar, double lo,
                                 double hi, unsigned mode, double r1, double r2)
    : RooAbsReal(name, title),
      xvar_("xvar", "Variable", this, xvar),
      vsmooth_par_(1.0),
      lo_(lo),
      hi_(hi),
      mode_(mode) {
        if (mode == 1 || mode == 2) {
          i0_.r_0_0_ = lo;
          i0_.r_0_1_ = -1. * lo;
          i0_.r_0_2_ = 0.;
          i0_.r_1_0_ = 0.;
          i0_.r_1_1_ = (hi - lo) / 2.;
          i0_.r_1_2_ = (hi + lo) * (15. / 8.);

          i1_.r_0_0_ = 0;
          i1_.r_0_1_ = (hi - lo) / 2.;
          i1_.r_0_2_ = (hi + lo) * (15. / 8.);;
          i1_.r_1_0_ = hi;
          i1_.r_1_1_ = hi;
          i1_.r_1_2_ = 0.;

          if (mode == 2) {
            i0_.r_1_1_ = r1;
            i0_.r_1_2_ = r2;
            i1_.r_0_1_ = r1;
            i1_.r_0_2_ = r2;
            mode_ = 1;
          }

          i0_.init();
          i1_.init();
        }
      }

RooCMSInterpVar::RooCMSInterpVar(const RooCMSInterpVar &other, const char *newname)
    : RooAbsReal(other, newname), xvar_("xvar", this, other.xvar_) {}

RooCMSInterpVar::~RooCMSInterpVar() {}

TObject *RooCMSInterpVar::clone(const char *newname) const {
  return new RooCMSInterpVar(*this, newname);
}

Double_t RooCMSInterpVar::evaluate() const {
  if (mode_ == 0) {
    return 0.5 * xvar_ * ((hi_ - lo_) + smoothStepFunc(xvar_) * (hi_ + lo_));
  }
  if (mode_ == 1 && xvar_ <= 0.) {
    double t = xvar_ + 1;
    return i0_.eval(t);
  }
  if (mode_ == 1 && xvar_ > 0.) {
    double t = xvar_;
    return i1_.eval(t);
  }
  return 0;
}

double RooCMSInterpVar::NewInterp::Bern(unsigned p, unsigned i, double t) const {
  return TMath::Binomial(p, i) * std::pow(t, i) * std::pow(1. - t, p - i);
}

void RooCMSInterpVar::NewInterp::init() {
  sig_ = 5.0;
  std::cout << "a) sig := " << sig_ << "\n";
  del_ = std::pow(r_0_1_ + r_1_1_, 2) - (r_1_0_ - r_0_0_) * (r_1_2_ - r_0_2_);
  std::cout << "b) i. del := " << del_ << "\n";

  if (del_ > 0.) {
    double alt = 1. + (r_0_1_ + r_1_1_ + std::sqrt(del_)) / (r_1_0_ - r_0_0_);
    std::cout << "   ii.  sig = max(sig, " << alt << ") = ";
    sig_ = std::max(sig_, alt);
    std::cout << sig_ << "\n";
  }

  if (r_0_1_ > 0.) {
    double alt = 1. - r_0_2_ / r_0_1_;
    std::cout << "   iii. sig = max(sig, " << alt << ") = ";
    sig_ = std::max(sig_, alt);
    std::cout << sig_ << "\n";
  }
  if (r_1_1_ > 0.) {
    double alt = 1. + r_1_2_ / r_1_1_;
    std::cout << "   iv.  sig = max(sig, " << alt << ") = ";
    sig_ = std::max(sig_, alt);
    std::cout << sig_ << "\n";
  }

  c0_ = r_0_0_;
  c1_ = r_0_0_ + (1. / sig_) * r_0_1_;
  c2_ = r_0_0_ + (2. / sig_) * r_0_1_ + (1. / (sig_ * (sig_ - 1.))) * r_0_2_;
  c3_ = r_1_0_ - (2. / sig_) * r_1_1_ + (1. / (sig_ * (sig_ - 1.))) * r_1_2_;
  c4_ = r_1_0_ - (1. / sig_) * r_1_1_;
  c5_ = r_1_0_;

  w0_ = 1.;
  w1_ = (sig_ - 1.) / 4.;
  w2_ = ((sig_ - 1.) * (sig_ - 2.)) / 12.;
  w3_ = (sig_ - 1.) / 4.;
  w4_ = 1.;

  wb0_ = 1.;
  wb1_ = sig_ / 5.;
  wb2_ = (sig_ * (sig_ - 1.)) / 20.;
  wb3_ = (sig_ * (sig_ - 1.)) / 20.;
  wb4_ = sig_ / 5.;
  wb5_ = 1.;
}

void RooCMSInterpVar::NewInterp::initConcave() {
  sig_ = 5.0;
  std::cout << "a) sig := " << sig_ << "\n";
  del_ = std::pow(r_0_1_ - r_1_1_ - r_0_2_ / 2., 2) - (r_1_0_ - r_0_0_ - r_0_1_) * (r_1_2_ + 2. * r_0_2_);
  std::cout << "b) i. del := " << del_ << "\n";

  if (del_ > 0.) {
    double alt = 1. + (-r_0_1_ + r_1_1_ + r_0_2_ / 2. + std::sqrt(del_)) / (r_1_0_ - r_0_0_ - r_0_1_);
    std::cout << "   ii.  sig = max(sig, " << alt << ") = ";
    sig_ = std::max(sig_, alt);
    std::cout << sig_ << "\n";
  }
  double del1 = std::pow(r_0_1_ - r_1_1_ - r_1_2_ / 2., 2) - (r_1_1_ - r_1_0_ - r_0_0_) * (r_0_2_ + 2. * r_1_2_);
  if (del1 > 0.) {
    double alt = 1. + (-r_0_1_ + r_1_1_ + r_1_2_ / 2. + std::sqrt(del1)) / (r_1_1_ - r_1_0_ - r_0_0_);
    std::cout << "   iii. sig = max(sig, " << alt << ") = ";
    sig_ = std::max(sig_, alt);
    std::cout << sig_ << "\n";
  }

  c0_ = r_0_0_;
  c1_ = r_0_0_ + (1. / sig_) * r_0_1_;
  c2_ = r_0_0_ + (2. / sig_) * r_0_1_ + (1. / (sig_ * (sig_ - 1.))) * r_0_2_;
  c3_ = r_1_0_ - (2. / sig_) * r_1_1_ + (1. / (sig_ * (sig_ - 1.))) * r_1_2_;
  c4_ = r_1_0_ - (1. / sig_) * r_1_1_;
  c5_ = r_1_0_;

  w0_ = 1.;
  w1_ = (sig_ - 1.) / 4.;
  w2_ = ((sig_ - 1.) * (sig_ - 2.)) / 12.;
  w3_ = (sig_ - 1.) / 4.;
  w4_ = 1.;

  wb0_ = 1.;
  wb1_ = sig_ / 5.;
  wb2_ = (sig_ * (sig_ - 1.)) / 20.;
  wb3_ = (sig_ * (sig_ - 1.)) / 20.;
  wb4_ = sig_ / 5.;
  wb5_ = 1.;
}

double RooCMSInterpVar::NewInterp::eval(double t) const {
  double P = wb0_ * c0_ * Bern(5, 0, t) +
             wb1_ * c1_ * Bern(5, 1, t) +
             wb2_ * c2_ * Bern(5, 2, t) +
             wb3_ * c3_ * Bern(5, 3, t) +
             wb4_ * c4_ * Bern(5, 4, t) +
             wb5_ * c5_ * Bern(5, 5, t);

  double Q = w0_ * Bern(4, 0, t) +
             w1_ * Bern(4, 1, t) +
             w2_ * Bern(4, 2, t) +
             w3_ * Bern(4, 3, t) +
             w4_ * Bern(4, 4, t);

  double R = P / Q;
  return R;
}


ClassImp(RooCMSInterpVar)
