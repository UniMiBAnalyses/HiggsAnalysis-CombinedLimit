/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

/**
 * @example MCMC.cpp
 * This is an example of how to perform a Bayesian Markov Chain Monte Carlo analysis with HEPfit and BAT.
 *
 */

#include <iostream>
#include <RooRealVar.h>
#include <RooMinimizer.h>
#include "TFile.h"
#include "TTree.h"
// #include <HEPfit.h>
// #include <BAT/BCParameter.h>
/* Necessary if MPI support is enabled during compilation. */
#include "HiggsAnalysis/CombinedLimit/interface/RooHEPFitWrapper.h"
int main(int argc, char** argv)
{

    RooHEPFitWrapper wrapper("HEPFitWrapper", "", "model.conf");
    wrapper.InitModel();
    auto vars = wrapper.InitVars();
    vars->Print("v");

    auto mHl = (RooRealVar*)vars->find("mHl");


    RooMinimizer minim(wrapper);
    minim.setVerbose(false);
    minim.setEps(0.1);
    minim.minimize("Minuit2", "migrad");

    TFile fout("hepfit_test.root", "RECREATE");
    TTree tout("limit", "limit");
    float POI = mHl->getVal();
    float deltaNLL = 0;
    float quantileExpected = 1;
    double nll0 = wrapper.getVal();

    tout.Branch("mHl", &POI, "mHl/f");
    tout.Branch("deltaNLL", &deltaNLL, "deltaNLL/f");
    tout.Branch("quantileExpected", &quantileExpected, "quantileExpected/f");

    tout.Fill();

    unsigned points = 80;
    double width = (mHl->getMax() - mHl->getMin()) / double(points);
    double r = mHl->getMin() + 0.5 * width;

    auto snap = vars->snapshot();

    for (unsigned ip = 0; ip < points; ++ip) {
        vars->assignValueOnly(*snap);
        mHl->setVal(r);
        mHl->setConstant(true);
        // minim.minimize("Minuit2", "migrad");
        POI = r;
        deltaNLL = wrapper.getVal() - nll0;
        std::cout << "var = " << r << "; nll0 = " << nll0 << ", deltaNLL = " << deltaNLL << "\n";
        tout.Fill();
        r += width;
    }

    tout.Write();
    fout.Close();
    // minim.save()->Print();
//     /* Necessary if MPI support is enabled during compilation. */
// #ifdef _MPI
//     MPI::Init();
//     int rank = MPI::COMM_WORLD.Get_rank();
//     MPI::Status status;
// #else
//     /* In our MPI implementation the process with rank 0 is the master. */
//     int rank = 0;
// #endif
    
//     try {
        
//         if(argc != 3){
//             /* Print usage and exit. */
//             if (rank == 0) std::cout << "\nusage: " << argv[0] << " ModelConf.conf --noMC/MonteCarlo.conf\n" << std::endl;
//             return EXIT_SUCCESS;
//         }

//         /* Define the model configuration file.                                */
//         /* Here it is passed as the first argument to the executable. The      */
//         /* model configuration file provides the values with errors for the    */
//         /* mandatory model parameters, as well as the list of observables,     */
//         /* observables2D, correlated Gaussian observables.                     */
//         /* See documentation for details.                                      */
//         std::string ModelConf = argv[1];
        
//         /* Define the Monte Carlo configuration file.                         */
//         /* Here it is passed as the second argument to the executable. The    */
//         /* Monte Carlo configuration file provides the parameters used in the */
//         /* Monte Carlo run. See documentation for details.                    */
//         std::string MCMCConf = argv[2];
        
//         /* Define the ROOT output file (w/o extension, empty string will set it to MCout) */
//         std::string FileOut = "";
        
//         /* Define the optional job tag. */
//         std::string JobTag = "";
        
//         /* Create objects of the classes ModelFactory and ThObsFactory */
//         ThObsFactory ThObsF;
//         ModelFactory ModelF;

//         /* register user-defined model named ModelName defined in class ModelClass using the following syntax: */
//         /* ModelF.addModelToFactory(ModelName, boost::factory<ModelClass*>() ) */
        
//         /* register user-defined ThObservable named ThObsName defined in class ThObsClass using the following syntax: */
//         /* ThObsF.addObsToFactory(ThObsName, boost::factory<ThObsClass*>() )*/
        
//         /* Create an object of the class MonteCarlo. */
//         MonteCarlo MC(ModelF, ThObsF, ModelConf, MCMCConf, FileOut, JobTag);
        
//          Do a test run if you wish to see the values of the observables      
//         /* and the correlated Gaussian observables defined in the model        */
//         /* configuration file computed with the central value of the mandatory */
//         /* parameters defined in the same file.                                */
//         if (MCMCConf.compare("--noMC") == 0) {

//             MC.TestRun(rank);
//             std::vector<double> newpars;
//             for (unsigned ip = 0; ip < MC.modelParameters().size(); ++ip) {
//                 ModelParameter const& param = MC.modelParameters()[ip];
//                 if (!param.IsFixed()) {
//                     std::cout << param.getname() << " = " << param.getave() << " (" << param.IsFixed() << ")\n";
//                     newpars.push_back(param.getave());
//                 }
//             }
//             std::cout << "GetNParameters: " << MC.getMCEngine()->GetNParameters() << "\n";
//             for (unsigned ip = 0; ip < MC.getMCEngine()->GetParameters().Size(); ++ip) {
//                 std::cout << MC.getMCEngine()->GetParameters()[ip]->GetName() << "\n";
//             }
//             newpars[3] = 98.0;
//             std::cout << "NLL = " << MC.getMCEngine()->Function_h(newpars) << "\t"
//                       << MC.getMCEngine()->LogLikelihood(newpars) << "\t"
//                       << MC.getMCEngine()->LogAPrioriProbability(newpars) << "\n";
//             newpars[3] = 128.0;
//             std::cout << "NLL = " << MC.getMCEngine()->Function_h(newpars) << "\t"
//                       << MC.getMCEngine()->LogLikelihood(newpars) << "\t"
//                       << MC.getMCEngine()->LogAPrioriProbability(newpars) << "\n";
//             for (unsigned ip = 0; ip < MC.modelParameters().size(); ++ip) {
//                 ModelParameter const& param = MC.modelParameters()[ip];
//                 if (!param.IsFixed()) {
//                     std::cout << param.getname() << " = " << param.getave() << " (" << param.IsFixed() << ")\n";
//                 }
//             }
//         }
//         /* Initiate the Mote Carlo run. */
//         else {
//             MC.Run(rank);
//         }

//     /* Necessary if MPI support is enabled during compilation. */
// #ifdef _MPI
//         MPI::Finalize();
// #endif
        
//         return EXIT_SUCCESS;
//     } catch (const std::runtime_error& e) {
//         std::cerr << e.what() << std::endl;
//         return EXIT_FAILURE;
//     }
}
