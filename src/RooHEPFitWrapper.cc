#include "RooFit.h"

#include "Riostream.h"
#include <math.h>

#include "../interface/RooHEPFitWrapper.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooMath.h"
#include "RooConstVar.h"
#include "RooListProxy.h"
#include "TArrayI.h"
#include <HEPfit.h>
#include <BAT/BCParameter.h>

using namespace std;

ClassImp(RooHEPFitWrapper);

//_____________________________________________________________________________
RooHEPFitWrapper::RooHEPFitWrapper(
    const char* name, const char* title, std::string const& conf)
    : RooAbsReal(name, title),
      _init(false),
      _conf(conf),
      _x("x", "Observables", this, kTRUE, kFALSE),
      _MC(nullptr) {
}

//_____________________________________________________________________________
RooHEPFitWrapper::RooHEPFitWrapper(const RooHEPFitWrapper& other,
                                       const char* name)
    : RooAbsReal(other, name),
      _init(false),
      _conf(other._conf),
      _x("x", this, other._x),
      _MC(nullptr) {}

//_____________________________________________________________________________
Double_t RooHEPFitWrapper::evaluate() const {

  std::vector<double> point(_x.getSize());
  for (int i = 0; i < _x.getSize(); ++i) {
    point[i] = ((RooAbsReal*)(_x.at(i)))->getVal();
  }
  double nll = -1 * _MC->getMCEngine()->Function_h(point);
  return nll;
}


void RooHEPFitWrapper::printMultiline(std::ostream& os, Int_t contents,
                                 Bool_t verbose, TString indent) const {
  RooAbsReal::printMultiline(os, contents, verbose, indent);
}

void RooHEPFitWrapper::InitModel() {
  if (_init) {
    std::cout << "Already initialised the model\n";
    return;
  }

  if (_MC  == nullptr) {
    /* Define the ROOT output file (w/o extension, empty string will set it to MCout) */
    std::string FileOut = "";
    /* Define the optional job tag. */
    std::string JobTag = "";
    /* Create objects of the classes ModelFactory and ThObsFactory */
    ThObsFactory ThObsF;
    ModelFactory ModelF;
    std::string MCMCConf = "MonteCarlo.conf";
    /* Create an object of the class MonteCarlo. */
    _MC = new MonteCarlo(ModelF, ThObsF, _conf, MCMCConf, FileOut, JobTag);
    _MC->TestRun(0);
  }

}

RooArgList * RooHEPFitWrapper::InitVars() {
  RooArgList * res = new RooArgList();
  for (unsigned ip = 0; ip < _MC->modelParameters().size(); ++ip) {
      ModelParameter const& param = _MC->modelParameters()[ip];
      if (!param.IsFixed()) {
          std::cout << param.getname() << " = " << param.getave() << " (" << param.getmin() << ", " << param.getmax() << ")\n";
          RooRealVar *var = new RooRealVar(TString(param.getname()), "", param.getave(), param.getmin(), param.getmax());
          res->addOwned(*var);
          _x.add(*var);
          var->Print();
      }
  }
  return res;
}

