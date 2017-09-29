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
#ifndef ROO_TAYLOREXPANSION
#define ROO_TAYLOREXPANSION

#include "RooAbsPdf.h"
#include "RooListProxy.h"
#include "TVectorD.h"

class RooRealVar;
// class RooFitResult;

#include <map>

class RooTaylorExpansion : public RooAbsReal {
 public:
  RooTaylorExpansion(){};
  RooTaylorExpansion(const char* name, const char* title,
                     const RooArgList& x, const RooArgList& x0, unsigned int order,
                     const TVectorD& terms);

  RooTaylorExpansion(const RooTaylorExpansion& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const {
    return new RooTaylorExpansion(*this, newname);
  }
  inline virtual ~RooTaylorExpansion() {}

  // Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const
  // char* rangeName=0) const ; Double_t analyticalIntegral(Int_t code, const
  // char* rangeName=0) const ;

  // Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars,
  // Bool_t staticInitOK=kTRUE) const; void initGenerator(Int_t code) ; void
  // generateEvent(Int_t code);


 protected:
  // void decodeCode(Int_t code, std::vector<int>& map1, std::vector<int>& map2)
  // const; AnaIntData& anaIntData(Int_t code) const ; GenData& genData(Int_t
  // code) const ;

  // mutable std::map<int, AnaIntData> _anaIntCache;  //!
  // mutable std::map<int, GenData> _genCache;        //!

  // mutable std::vector<BitBlock> _aicMap;  //!

  RooListProxy _x;
  RooListProxy _x0;
  unsigned int _order;
  TVectorD _terms;

  Double_t evaluate() const;

 private:
  ClassDef(RooTaylorExpansion,
           1)  // Multivariate Gaussian PDF with correlations
};

#endif
