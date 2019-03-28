#ifndef HiggsAnalysis_CombinedLimit_RooCMSInterpVar_h
#define HiggsAnalysis_CombinedLimit_RooCMSInterpVar_h

#include <RooAbsReal.h>
#include <RooRealProxy.h>

class RooCMSInterpVar : public RooAbsReal {
 public:
  RooCMSInterpVar() {}
  RooCMSInterpVar(const char *name, const char *title, RooAbsReal &xvar, double lo, double hi, unsigned mode, double r1, double r2);
  RooCMSInterpVar(const RooCMSInterpVar &other, const char *newname);
  ~RooCMSInterpVar();

  TObject *clone(const char *newname) const;

 protected:
  Double_t evaluate() const;

  inline double smoothStepFunc(double x) const {
    if (fabs(x) >= vsmooth_par_) return x > 0 ? +1 : -1;
    double xnorm = x / vsmooth_par_;
    double xnorm2 = xnorm * xnorm;
    return 0.125 * xnorm * (xnorm2 * (3. * xnorm2 - 10.) + 15);
  }

 private:

  struct NewInterp {
    double sig_;
    double del_;
    double r_0_0_;
    double r_0_1_;
    double r_0_2_;
    double r_1_0_;
    double r_1_1_;
    double r_1_2_;
    double c0_;
    double c1_;
    double c2_;
    double c3_;
    double c4_;
    double c5_;
    double w0_;
    double w1_;
    double w2_;
    double w3_;
    double w4_;
    double wb0_;
    double wb1_;
    double wb2_;
    double wb3_;
    double wb4_;
    double wb5_;
    void init();
    void initConcave();
    double eval(double t) const;
    double Bern(unsigned p, unsigned i, double t) const;
  };

  RooRealProxy xvar_;
  double vsmooth_par_;
  double lo_;
  double hi_;
  unsigned mode_;

  NewInterp i0_;
  NewInterp i1_;



  ClassDef(RooCMSInterpVar, 1)
};

#endif
