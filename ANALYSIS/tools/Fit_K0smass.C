#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooArgSet.h"

#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "My_double_CB/My_double_CB.cxx"

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
    {"2016preVFP_Psi", "../PRELIMINARY/outRoot/RecoDecay_Psi2S_UL16preVFP.root"},
    {"2016_Psi", "../PRELIMINARY/outRoot/RecoDecay_Psi2S_UL16.root"},
    {"2017_Psi", "../PRELIMINARY/outRoot/RecoDecay_Psi2S_UL17.root"},
    {"2018_Psi", "../PRELIMINARY/outRoot/RecoDecay_Psi2S_UL18.root"},
    {"run2_Psi", "../PRELIMINARY/outRoot/RecoDecay_Psi2S_ULrun2.root"}


};


void AddModel(RooWorkspace*, float, float);
void AddData(RooWorkspace* , string, float , float);
void FitToData(RooWorkspace* , string, float , float);
void SaveResults(RooWorkspace* , string, float , float);


void Fit_K0smass(string dataset = "2017_X"){

    // add the model
    Float_t lowMB  = 0.4; 
    Float_t highMB = 0.6;     
    RooWorkspace* wspace = new RooWorkspace("myWS");

    AddModel(wspace, lowMB, highMB);
    AddData(wspace, dataset, lowMB, highMB);

    wspace->Print();

    FitToData(wspace, dataset, lowMB, highMB);
    SaveResults(wspace, dataset, lowMB, highMB);    


}


void AddModel(RooWorkspace* ws, float low, float high){
   
    // variable for the observable 
    RooRealVar M_K0s("M_K0s", "K0s reco mass", low, high, "GeV");
    cout << " ... signal model for M(K0s)" << endl;

    // Double Gaussian 
    RooRealVar K0smass ("K0smass" , "", 0.495, 0.493, 0.499); 
    RooRealVar K0s_sigma1("K0s_sigma1", "", 0.005, 0.0001  , 0.01);
    RooRealVar K0s_sigma2("K0s_sigma2", "", 0.01, 0.005  , 0.05);
    
    RooGaussian G1signal("G1signal", "", M_K0s, K0smass, K0s_sigma1);
    RooGaussian G2signal("G2signal", "", M_K0s, K0smass, K0s_sigma2);
    // Double Side-CrystalBall
    RooRealVar K0s_sigma("K0s_sigma", "", 0.005, 0.0005  , 0.01);
    RooRealVar alphaL("alphaL", "", 1., 0., 10);
    RooRealVar nL("nL", "", 10, -10, 10);
    RooRealVar alphaR("alphaR", "", 1., 0., 10);
    RooRealVar nR("nR", "", 10, -10, 20);
    
    My_double_CB SGNmodel("SGNmodel", "", M_K0s, K0smass, K0s_sigma, alphaL, nL, alphaR, nR);

    RooRealVar f("f", "", 0., 1.);
    RooRealVar yield("yield", "", 0., 10000);
    //RooAddPdf SGNmodel("SGNmodel", "SGNmodel", RooArgList(G1signal, G2signal), f);
    //RooAddPdf SGNmodel("SGNmodel", "SGNmodel", CB2signal);


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
    RooRealVar *M_K0s = ws->var("M_K0s");
    RooDataSet data("data", "data", *M_K0s, Import(*tree));
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
    RooRealVar *M_K0s = ws->var("M_K0s");
    
    RooFitResult* results = model->fitTo(*data, Save());
    

    // create the plot
    int Nbins = 40;
    RooPlot *frame = M_K0s->frame(Title("MC-matching reco K0s candidates " + TString(dataset)));
    data->plotOn(frame, Binning(Nbins));
    model->plotOn(frame, LineColor(kGreen));
    double Chi2 = frame->chiSquare(6);
    data->plotOn(frame, Binning(Nbins));

    // save fit results on file
    TString outFileParams = "./outRoot/SGN_BKGregion/K0smassFit_" +dataset+".root"; 
    TFile* outFile = new TFile(outFileParams, "RECREATE");
    frame->Write("FitPlot");
    results->Write("FitParams");
    outFile->Close();
    cout << " ... fit results saved in " << outFileParams << endl;

    // compute and add to the wspace Chi2
    cout << " === Chi2 " << Chi2 << endl;
    RooRealVar FitChi2("Chi2", "SGNmodel-fit chi2", Chi2); 
    ws->import(FitChi2);
    
    RooCurve* FitCurve = (RooCurve*)frame->getObject(1); 

    TCanvas* cdata = new TCanvas("cdata","data fit", 800, 600);
    TLegend* legend = new TLegend(0.60, 0.75, 0.89, 0.89);

    gPad->SetMargin(0.12,0.1,0.12,0.1);
    frame->GetXaxis()->SetTitle("M(K^{0}_{s})");
    frame->Draw();

    legend->SetBorderSize(0);
    legend->AddEntry("data", "MC-matching cand.", "P");
    legend->AddEntry(FitCurve, "TOTAL fit", "LP");
    legend->Draw();

    // save the outputs
    TString out_name = "/eos/user/c/cbasile/www/B0toX3872K0s/RECO_LEVEL/MC_MassDistFit/Fit_K0smass_" + TString(dataset); 
    cdata->SaveAs(out_name + ".png");
    cdata->SaveAs(out_name + ".pdf");

    cdata->SetLogy(1);
    cdata->SaveAs(out_name + "_logY.png");
    cdata->SaveAs(out_name + "_logY.pdf");
    

} //fitToData()



void SaveResults(RooWorkspace* ws, string dataset, float Mlow, float Mhigh){

    
    RooAbsPdf* model      = ws->pdf("SGNmodel");
    RooDataSet* data = (RooDataSet*) ws->data("data");

    RooRealVar* M_K0s     = ws->var("M_K0s");
    RooRealVar* K0smass   = ws->var("K0smass");
    RooRealVar* K0s_sigma = ws->var("K0s_sigma");

    RooRealVar* chi2      = ws->var("Chi2");
    RooRealVar* data_N= ws->var("data_N");

	// SIDEBANDS
    double norm     = model->createIntegral(*M_K0s, NormSet(*M_K0s))->getVal();
    double Km       = K0smass->getVal();
	float SIGMA       = K0s_sigma->getVal();
	float SIGMA_ERR   = K0s_sigma->getError();
    cout << " === sigma_tot = " << SIGMA << " +/- " << SIGMA_ERR << endl;
    
    int SR_Ns= 3;
    float SGNeff = data->sumEntries(Form("M_K0s>%f&&M_K0s<%f", Km - SR_Ns*SIGMA, Km + SR_Ns*SIGMA))/ data_N->getVal();
    cout << " === signal-efficiency " << SGNeff << endl;
    
	string outFileParams = "./outRoot/SGN_BKGregion/K0smassFit_" +dataset+".txt";
	ofstream outFile;
	outFile.open(outFileParams);
	if(!outFile.is_open()) std::cout << "ERROR cannot open file" << outFileParams << std::endl;
	outFile << "### K0s signal region" << std::endl;	
    
	outFile << "sigma_tot"          << "\t" << SIGMA << " +/- " << SIGMA_ERR << "\n";
	outFile << "SR_low"             << "\t" << Km - SR_Ns*SIGMA << "\n"; 
	outFile << "SR_high"            << "\t" << Km + SR_Ns*SIGMA << "\n"; 
	outFile << "\n#########\n";
	outFile << "chi-square"            << "\t" << chi2->getVal() << "\n"; 
	outFile << "signal efficiency \t"<< SGNeff << "\n"; 
    
	outFile.close();
    cout << " ... signal region saved in " << outFileParams << endl;


}//SaveResults()

