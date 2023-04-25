#ifndef OptimizerMVA_h
#define OptimizerMVA_h

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "tools/RootPlots.h"


// RooFit
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooFFTConvPdf.h"
#include "RooGenericPdf.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooCurve.h"
#include "RooWorkspace.h"
using namespace RooFit;



#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
 
#include "TFile.h"
#include "TTree.h"
#include "TText.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStopwatch.h"
 
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include <TStyle.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

class MVAoptimizer{

public:

    // constructor & destructor
    MVAoptimizer(const TString& input, const float& BDTcut,const float& Mcut,const int& year, const TString & channel = "X_3872");
    ~MVAoptimizer();

    //getters
    float get_BDT_cut(){return BDToutCut_;}
    float get_Mpipi_cut(){return MpipiCut_;}

    // setters
    void set_SBregions();
    void set_Blind(const bool& blindness){is_blind_ = blindness;}
    void set_LumiNorm(const float &factor){LumiNorm_ = factor;}
    void set_BDT_cut(float cut){BDToutCut_ = cut;}
    void set_Mpipi_cut(float cut){MpipiCut_ = cut;}
    void set_selection();

    // optimization tools
    double Nbkg_extraction();
    double Nsgn_extraction();
    double EFFsgn_extraction();
    double PunziSign(double * PSerr);

    int  makeSGNvsBKGplot();
    void Plot2D_BDToutMpipi();

    void CMSxxx(TCanvas*);

private:
    // dataset
    int year_;
    float LumiNorm_;
    TString channel_;
    bool is_blind_;

    // input data
    TString inPath_;
    TFile*  inFile_;
    TTree*  inTree_;

    // output path
    TString outPath_;

    // plot interval
    double MB_low, MB_high;
    double MJpsiPiPi_low, MJpsiPiPi_high;
    double MK0s_low, MK0s_high;
    double MPiPi_low, MPiPi_high;
    double BDTout_low, BDTout_high;

    // cut thresholds
    float BDToutCut_, MpipiCut_;

    // signal and background edges
	double B0_SR_low, B0_SR_high;
	double JpsiPiPi_SR_low, JpsiPiPi_SR_high;
	double K0s_SR_low, K0s_SR_high;
	double B0_lSB_low, B0_lSB_high, B0_rSB_low, B0_rSB_high;

    // selection
    TString SGNselection, DATAselection, SRselection;
    TString selection_;

    // RooFit workspace
    RooWorkspace * wsFit = new RooWorkspace("wsFit", "workspace");


}; // class MVAoptimizer
#endif
