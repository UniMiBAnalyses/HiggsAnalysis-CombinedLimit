/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_HEPFITWRAPPER
#define ROO_HEPFITWRAPPER

#include <string>
#include "RooAbsPdf.h"
#include "RooListProxy.h"
// #include <HEPfit.h>

class RooRealVar;
class MonteCarlo;

#include <map>

class RooHEPFitWrapper : public RooAbsReal {
 public:
  RooHEPFitWrapper(){};
  RooHEPFitWrapper(const char* name, const char* title, std::string const& conf);

  RooHEPFitWrapper(const RooHEPFitWrapper& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const {
    return new RooHEPFitWrapper(*this, newname);
  }
  inline virtual ~RooHEPFitWrapper() {}

  void printMultiline(std::ostream& os, Int_t contents,
                                   Bool_t verbose, TString indent) const;

  void InitModel();
  RooArgList * InitVars();

 protected:
  bool _init;
  std::string _conf;
  RooListProxy _x;

  MonteCarlo * _MC; //! not to be serialized

  Double_t evaluate() const;



 private:
  ClassDef(RooHEPFitWrapper,
           1)  // Multivariate Gaussian PDF with correlations
};

#endif
