#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooArgSet.h"

#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooCBShape.h"

#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"

#include "RooStats/SPlot.h"

#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"

#include "TSystem.h"

#include <map>
#include <iostream>


// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using namespace std;


map<string, TString> MC_rootFile{

    {"2016preVFP_X", "../PRELIMINARY/outRoot/RecoDecay_X3872_UL16preVFP.root"},
    {"2016_X", "../PRELIMINARY/outRoot/RecoDecay_X3872_UL16.root"},
    {"2017_X", "../PRELIMINARY/outRoot/RecoDecay_X3872_UL17.root"},
    {"2018_X", "../PRELIMINARY/outRoot/RecoDecay_X3872_UL18.root"},
    {"run2_X", "../PRELIMINARY/outRoot/RecoDecay_X3872_Run2.root"},
    {"2016preVFP_Psi", "../PRELIMINARY/outRoot/RecoDecay_Psi2S_UL16preVFP.root"},
    {"2016_Psi", "../PRELIMINARY/outRoot/RecoDecay_Psi2S_UL16.root"},
    {"2017_Psi", "../PRELIMINARY/outRoot/RecoDecay_Psi2S_UL17.root"},
    {"2018_Psi", "../PRELIMINARY/outRoot/RecoDecay_Psi2S_UL18.root"},
    {"run2_Psi", "../PRELIMINARY/outRoot/RecoDecay_Psi2S_Run2.root"}


};


void AddModel(RooWorkspace*, float, float);
void AddData(RooWorkspace* , string, float , float);
void FitToData(RooWorkspace* , string, float , float);
void SaveResults(RooWorkspace* , string, float , float);


void Fit_B0mass(string dataset = "2017_X"){

    // add the model
    Float_t lowMB  = 5.0; 
    Float_t highMB = 5.6;     
    RooWorkspace* wspace = new RooWorkspace("myWS");

    AddModel(wspace, lowMB, highMB);
    AddData(wspace, dataset, lowMB, highMB);

    wspace->Print();

    FitToData(wspace, dataset, lowMB, highMB);
    SaveResults(wspace, dataset, lowMB, highMB);    


}


void AddModel(RooWorkspace* ws, float low, float high){
   
    // variable for the observable 
    RooRealVar M_B0("M_B0", "B0 reco mass", low, high, "GeV");
    cout << " ... signal model for M(B0)" << endl;

    // Crystal Ball
    RooRealVar B0mass ("B0mass" , "", 5.279, 5.275, 5.285);
    RooRealVar B0width("B0width", "", 0.01, 0.005  , 0.1);
    RooRealVar alpha("alpha", "", 2., 0.,  20.);
    RooRealVar     N(    "N", "", 1., 0.0001,  100.);
    
    RooCBShape CBsignal("CBsignal", "", M_B0, B0mass, B0width, alpha, N);
    RooProdPdf SGNpdfCB("SGNpdfCB", "SGNpdfCB", CBsignal ); 
    // Gaussian
    RooRealVar B0sigma("B0sigma", "", 0.02, 0.0  , 2.);
    
    RooGaussian Gsignal("Gsignal", "", M_B0, B0mass, B0sigma);
    
    RooRealVar f("f", "", 0., 1.);
    RooAddPdf SGNmodel("SGNmodel", "SGNmodel", RooArgList(CBsignal, Gsignal), f);


    // import model in the workspace
    cout << " --> import the full model" << endl;
    ws->import(SGNmodel);

} //AddModel()



void AddData(RooWorkspace* ws, string dataset, float Mlow, float Mhigh){


	// load Tree-branch 
	TString tree_name = "RecoDecay_MCmatch";
	TFile* input_file = new TFile(MC_rootFile[dataset]);
	TTree* tree = (TTree*)input_file->Get(tree_name);
    if ( !tree ){
        std::cout<< "null pointer for TTree named " << tree_name << std::endl;
        exit(-1);
    }

    cout << " + got data from tree " << tree_name << " from file "<< MC_rootFile[dataset] << endl;
    RooRealVar *M_B0 = ws->var("M_B0");
    RooDataSet data("data", "data", *M_B0, Import(*tree));
    RooRealVar data_N("data_N", "", tree->GetEntries());
    data.Print();

    ws->import(data);
    ws->import(data_N);
    
    input_file->Close();


}//AddData()



void FitToData(RooWorkspace* ws, string dataset, float Mlow, float Mhigh){

    cout << " ... start fitting data ..." << std::endl;
    
    RooAbsPdf* model = ws->pdf("SGNmodel");
    RooDataSet* data = (RooDataSet*) ws->data("data");
    
    RooFitResult* results = model->fitTo(*data, Save());
    
    // save parameters values after the fit
    RooRealVar* M_B0     = ws->var("M_B0");
    RooRealVar* B0mass  = ws->var("B0mass");
    RooRealVar* B0sigma = ws->var("B0sigma");
    RooRealVar* B0width = ws->var("B0width");
    RooRealVar* alpha   = ws->var("alpha");
    RooRealVar* N       = ws->var("N");

    // create the plot
    RooPlot *frame = M_B0->frame(Title("MC-matching reco B0 candidates " + TString(dataset)));
    data->plotOn(frame, Binning(60));
    model->plotOn(frame, Components("Gsignal"), LineColor(kOrange), LineStyle(kDashed));
    model->plotOn(frame, Components("CBsignal"), LineColor(kRed), LineStyle(kDashed));
    model->plotOn(frame, LineColor(kAzure +1));
    double Chi2 = frame->chiSquare(5);
    data->plotOn(frame, Binning(60));

    // save fit results on file
    TString outFileParams = "./outRoot/SGN_BKGregion/B0massFit_" +dataset+".root"; 
    TFile* outFile = new TFile(outFileParams, "RECREATE");
    frame->Write("FitPlot");
    results->Write("FitParams");
    outFile->Close();
    cout << " ... fit results saved in " << outFileParams << endl;

    // compute and add to the wspace Chi2
    cout << " === Chi2 " << Chi2 << endl;
    RooRealVar FitChi2("Chi2", "SGNmodel-fit chi2", Chi2); 
    ws->import(FitChi2);
    
    RooCurve* FitG     = (RooCurve*)frame->getObject(1); 
    RooCurve* FitCB    = (RooCurve*)frame->getObject(2); 
    RooCurve* FitCurve = (RooCurve*)frame->getObject(3); 

    TCanvas* cdata = new TCanvas("cdata","data fit", 800, 600);
    TLegend* legend = new TLegend(0.60, 0.70, 0.89, 0.89);

    gPad->SetMargin(0.12,0.1,0.12,0.1);
    frame->GetXaxis()->SetTitle("M(B^{0})");
    frame->Draw();

    legend->SetBorderSize(0);
    legend->AddEntry("data", "MC-matching cand.", "P");
    legend->AddEntry(FitCurve, "TOTAL fit", "LP");
    legend->AddEntry(FitG, "Gaussian", "LP");
    legend->AddEntry(FitCB, "Crystal Ball", "LP");
    legend->Draw();

    // save the outputs
    TString out_name = "/eos/user/c/cbasile/www/B0toX3872K0s/RECO_LEVEL/MC_MassDistFit/Fit_B0mass_" + TString(dataset); 
    cdata->SaveAs(out_name + ".png");
    cdata->SaveAs(out_name + ".pdf");

    cdata->SetLogy(1);
    cdata->SaveAs(out_name + "_logY.png");
    cdata->SaveAs(out_name + "_logY.pdf");
    

} //FitToData()



void SaveResults(RooWorkspace* ws, string dataset, float Mlow, float Mhigh){

    
    RooAbsPdf* model    = ws->pdf("SGNmodel");
    RooDataSet* data    = (RooDataSet*) ws->data("data");

    RooRealVar* M_B0    = ws->var("M_B0");
    RooRealVar* B0mass  = ws->var("B0mass");
    RooRealVar* B0sigma = ws->var("B0sigma");
    RooRealVar* B0width = ws->var("B0width");
    RooRealVar* alpha   = ws->var("alpha");
    RooRealVar* N       = ws->var("N");
    RooRealVar* f       = ws->var("f");

    RooRealVar* chi2    = ws->var("Chi2");
    RooRealVar* data_N  = ws->var("data_N");

	// SIDEBANDS
    double norm     = model->createIntegral(*M_B0, NormSet(*M_B0))->getVal();
    double Bm       = B0mass->getVal();
    double wG       = 1-f->getVal();
    double sG       = B0sigma->getVal();
    double sG_err   = B0sigma->getError();
    double wCB      = f->getVal();
    double sCB      = B0width->getVal();
    double sCB_err  = B0width->getError();
	float SIGMA     = 2 * sqrt( (sG*sG*wG*wG  )  + (sCB*sCB*wCB*wCB));
	float SIGMA_ERR = 4./SIGMA * sqrt( (wG*wG*sG*sG_err)*(wG*wG*sG*sG_err) + (wCB*wCB*sCB*sCB_err)*(wCB*wCB*sCB*sCB_err));
    cout << " === sigma_tot = " << SIGMA << " +/- " << SIGMA_ERR << endl;
    
    int SB_Ns= 4;
    int SR_Ns= 3;
    float SGNeff  = data->sumEntries(Form("M_B0>%f&&M_B0<%f", Bm - SR_Ns*SIGMA, Bm + SR_Ns*SIGMA))/ data_N->getVal();
    float SGNcont = data->sumEntries(Form("(M_B0>%f&&M_B0<%f)||(M_B0>%f&&M_B0<%f)", Bm - 2*SB_Ns*SIGMA, Bm - SB_Ns*SIGMA, Bm + SB_Ns*SIGMA, Bm + 2*SB_Ns*SIGMA))/ data_N->getVal();
    cout << " === signal-efficiency " << SGNeff << endl;
    cout << " === signal-contamination " << SGNcont << endl;
    
	string outFileParams = "./outRoot/SGN_BKGregion/B0massFit_" +dataset+".txt";
	ofstream outFile;
	outFile.open(outFileParams);
	if(!outFile.is_open()) std::cout << "ERROR cannot open file" << outFileParams << std::endl;
	outFile << "### B0 sidebands and signal region" << std::endl;	
    
	outFile << "sigma_tot"          << "\t" << SIGMA << " +/- " << SIGMA_ERR << "\n";
	outFile << "SBl_low"            << "\t" << Bm - 2*SB_Ns*SIGMA << "\n";
	outFile << "SBl_high"           << "\t" << Bm -   SB_Ns*SIGMA << "\n";
	outFile << "SBr_low"            << "\t" << Bm +   SB_Ns*SIGMA << "\n";
	outFile << "SBr_high"           << "\t" << Bm + 2*SB_Ns*SIGMA << "\n\n";
	outFile << "SR_low"             << "\t" << Bm - SR_Ns*SIGMA << "\n"; 
	outFile << "SR_high"            << "\t" << Bm + SR_Ns*SIGMA << "\n"; 
	outFile << "\n#########\n";
	outFile << "chi-square"            << "\t" << chi2->getVal() << "\n"; 
	outFile << "signal efficiency \t"<< SGNeff << "\n"; 
	outFile << "signal contamination \t"<< SGNcont<< "\n"; 
    
	outFile.close();
    cout << " ... signal and sidebands region saved in " << outFileParams << endl;


}//SaveResults()

