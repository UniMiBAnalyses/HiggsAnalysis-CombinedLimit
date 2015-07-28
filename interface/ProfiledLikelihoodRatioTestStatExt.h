#ifndef HiggsAnalysis_CombinedLimit_ProfiledLikelihoodRatioTestStatExt_h
#define HiggsAnalysis_CombinedLimit_ProfiledLikelihoodRatioTestStatExt_h

#include <memory>
#include <vector>

class RooMinimizerOpt;
#include <RooAbsPdf.h>
#include <RooAbsData.h>
#include <RooArgSet.h>
#include <RooStats/TestStatistic.h>
#include "../interface/RooSimultaneousOpt.h"
#include "../interface/CachingNLL.h"

namespace nllutils {
    bool robustMinimize(RooAbsReal &nll, RooMinimizerOpt &minimizer, int verbosity=0, bool zeroPoint=false);
}

class ProfiledLikelihoodRatioTestStatOpt : public RooStats::TestStatistic {
    public:
        ProfiledLikelihoodRatioTestStatOpt(const RooArgSet &obs, RooAbsPdf &pdfNull, RooAbsPdf &pdfAlt, 
                const RooArgSet *nuisances, const RooArgSet & paramsNull = RooArgSet(), const RooArgSet & paramsAlt = RooArgSet(),
                int verbosity=0) ;
 
        virtual Double_t Evaluate(RooAbsData& data, RooArgSet& nullPOI) ;

        virtual const TString GetVarName() const {
            return TString::Format("-log(%s/%s)", pdfNull_->GetName(), pdfAlt_->GetName()); 
        }

        // Verbosity (default: 0)
        void setPrintLevel(Int_t level) { verbosity_ = level; }

        void SetupOutput(bool save_tree, bool save_data, double data_min_q,
                         double data_max_q);

        virtual ~ProfiledLikelihoodRatioTestStatOpt();

    private:
        RooAbsPdf *pdfNull_, *pdfAlt_;
        RooArgSet snapNull_, snapAlt_; 
        RooArgSet nuisances_; 
        std::auto_ptr<RooArgSet> paramsNull_, paramsAlt_;
        std::auto_ptr<RooAbsReal> nllNull_, nllAlt_;
        Int_t verbosity_;
        TFile *f_out_;
        TTree *t_out_;
        double save_data_min_;
        double save_data_max_;
        bool save_data_;
        bool save_tree_;
        double t_deltaNLL_;
        bool t_fit_nul_;
        bool t_fit_alt_;
        int t_index_;


        // create NLL. if returns true, it can be kept, if false it should be deleted at the end of Evaluate
        bool createNLL(RooAbsPdf &pdf, RooAbsData &data, std::auto_ptr<RooAbsReal> &nll) ;

        double minNLL(std::auto_ptr<RooAbsReal> &nll) ;
}; // TestSimpleStatistics


class ProfiledLikelihoodTestStatOpt : public RooStats::TestStatistic {
    public:
        enum OneSidedness { twoSidedDef = 0, oneSidedDef = 1, signFlipDef = 2 };

        ProfiledLikelihoodTestStatOpt(const RooArgSet & observables,
                RooAbsPdf &pdf, 
                const RooArgSet *nuisances, 
                const RooArgSet & params, const RooArgSet & poi, const RooArgList &gobsParams, const RooArgList &gobs, int verbosity=0, OneSidedness oneSided = oneSidedDef) ; 

        virtual Double_t Evaluate(RooAbsData& data, RooArgSet& nullPOI) ;
        virtual std::vector<Double_t> Evaluate(RooAbsData& data, RooArgSet& nullPOI, const std::vector<Double_t> &rVals) ;

        virtual const TString GetVarName() const { return "- log (#lambda)"; }

        // Verbosity (default: 0)
        void setPrintLevel(Int_t level) { verbosity_ = level; }

        void SetOneSided(OneSidedness oneSided) { oneSided_ = oneSided; }
    private:

        RooAbsPdf *pdf_;
        RooArgSet snap_, poi_; // snapshot of parameters, and of the subset which are POI
        std::auto_ptr<RooArgSet>  params_;
        RooArgSet                 nuisances_; // subset of params which are nuisances (not a snapshot)
        RooArgSet                 poiParams_; // subset of params which are POI (not a snapshot)
        std::auto_ptr<RooAbsReal> nll_;
        RooArgList gobsParams_, gobs_;
        Int_t verbosity_;
        OneSidedness oneSided_;

        // create NLL. if returns true, it can be kept, if false it should be deleted at the end of Evaluate
        bool createNLL(RooAbsPdf &pdf, RooAbsData &data) ;
        double minNLL(bool constrained, RooRealVar *r=0) ;
}; // TestSimpleStatistics


#endif
