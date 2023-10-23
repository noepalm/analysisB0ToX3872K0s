#include "../include/MVAoptimizer.h"
#include "../include/LumiConstants.h"
#include "RooCrystalBall.h"
using namespace LumiConstants; 

MVAoptimizer::MVAoptimizer(const TString& input, const float& BDTcut,const float& Mcut,const TString& year, const TString & channel){

    year_ = year;
    channel_ = channel;
	LumiNorm_ = McNormPerYear_Psi2S[year_]; 
	std::cout<< Form(" Analyze "+ year_ + " dataset for " + channel_ + " channel -> LumiNorm %.4f", LumiNorm_)<<std::endl;

    BDToutCut_ = BDTcut;
    MpipiCut_  = Mcut;
	is_blind_ = true;

    inPath_ = input;
	inTree_ = new TChain();

	load_files();
    set_SBregions();
    set_selection();

    if(channel_ == "Psi2S") outPath_ = "/eos/home-n/npalmeri/www/Analysis/Fit/"+ year_ + "/";
	else outPath_ = "/eos/home-n/npalmeri/www/Analysis/Fit/"+ year_ + "_test2/";
    //outPath_ = Form("/eos/user/c/cbasile/www/B0toX3872K0s/HLTchecks/");

}// MVAoptimizer()

MVAoptimizer::~MVAoptimizer(){
    inFile_->Close();
	wsFit->Delete();
}//~MVAoptimizer()


void MVAoptimizer::load_files(){

	std::array<TString, 5> eras = {"D", "E", "F", "G", ""}; //last for MC

	/*** WRITE TTREES WITH SELECTION TO FILE ***/
	TString outpath = "/eos/home-n/npalmeri/B0toX3872K0s/MVAresults/forFit/";

	// Add DATA trees to chain
	for(auto era : eras){

		// get right file name [use provided filename for path]
		TString file_name = inPath_;
		file_name.ReplaceAll(year_ + "D", year_ + era); 
		
    	if(!gSystem->AccessPathName(file_name)){

			// retrieve TTree
			TFile* f = TFile::Open(file_name);
			TTree* t = (TTree*)f->Get("CVtraining_" + year_ + era);

			// create output file to save TTree with selection
			TString outname = era == "" ? "MC" : "data";
			TFile* outf = TFile::Open(outpath + Form("%s%s_%s_%s.root", year_.Data(), era.Data(), outname.Data(), channel_.Data()), "RECREATE");

			// apply selection
			TString selection = era == "" ? "is_signal == 1" : "is_signal == 0";
			TTree* t_sel = t->CopyTree(selection);
			std::cout << " [+] wrote TTree CVtraining_" + year_ + era + " from " << file_name << " with " << t_sel->GetEntries() << " entries" << std::endl;

			// write and exit file
			t_sel->Write();
			outf->Close();

    	} else std::cout << " ERROR : cannot open input file " << file_name << std::endl;
	}

	/*** ADD FILES TO CHAIN ***/
	for(auto era : eras){
		TString outname = era == "" ? "MC" : "data";
		inTree_->Add(outpath + Form("%s%s_%s_%s.root/CVtraining_%s%s", year_.Data(), era.Data(), outname.Data(), channel_.Data(), year_.Data(), era.Data()));
	}

}

// // to be ALWAYS run after Nbkg_extraction(), Nsgn_extraction()
// void MVAoptimizer::create_pseudodataset(){

// 	// TODO: fix using 	https://root-forum.cern.ch/t/roofit-adding-two-datasets-one-scalled/5784 code

// 	// === ROOFIT SET UP === //
// 	// retrieve all vars
// 	RooRealVar* M_B0 = (RooRealVar*) wsFit->var("M_B0");
// 	RooRealVar* M_X3872 = (RooRealVar*) wsFit->var("M_X3872");
// 	RooRealVar* M_K0s = (RooRealVar*) wsFit->var("M_K0s");
// 	RooRealVar* M_PiPi = (RooRealVar*) wsFit->var("M_PiPi");
// 	RooRealVar* BDTout = (RooRealVar*) wsFit->var("BDTout");
// 	RooRealVar* is_signal = new RooRealVar("is_signal", "Flag for MC or data", 0, 1, ""); //needed for selection

// 	// === DATASET CREATION === //

// 	/*** BKG DATA, GENERATED ***/

// 	// retrieve BKGmodel from workspace
// 	RooAbsPdf* BKGmodel = wsFit->pdf("BKGmodel"); //NOTE: after full model fit, BKGmodel params are updated to latest values

// 	// retrieve n_bkg from full fit yield
// 	RooFitResult* ResBKGmodel = (RooFitResult*)wsFit->obj("fitresult_BKGmodel_data");
// 	double n_bkg_comb = ((RooRealVar*)(ResBKGmodel->floatParsFinal().find("n_comb")))->getVal();
// 	double n_bkg_jpsix = ((RooRealVar*)(ResBKGmodel->floatParsFinal().find("n_jpsix")))->getVal();
// 	double n_bkg = n_bkg_comb + n_bkg_jpsix;

// 	// generate dataset from BKGmodel
// 	RooDataSet* gen_data = BKGmodel->generate(*M_B0, NumEvents(n_bkg), Name("gen_data"));

// 	// retain only events IN BETWEEN sidebands
// 	gen_data = (RooDataSet*)gen_data->reduce(*M_B0, Form("M_B0 > %f && M_B0 < %f", B0_lSB_high, B0_rSB_low));

// 	/*** BKG DATA, SIDEBANDS ***/
// 	// [only select data with Cut(DATAselection), only select lSB, rSB range with CutRange("B0_lSB,B0_rSB")]
// 	TString DATA_selection = DATAselection;
// 	DATA_selection.Append(Form(" && M_X3872 > %f && M_X3872 < %f", JpsiPiPi_SR_low, JpsiPiPi_SR_high)); //only select signal events, i.e. tight around resonance
// 	DATA_selection.Append(Form(" && (M_B0 < %f || M_B0 > %f)", B0_lSB_high, B0_rSB_low));

// 	RooDataSet* sidebands_data = new RooDataSet("sidebands_data", "sidebands_data", RooArgSet(*M_B0, *M_X3872, *M_K0s, *M_PiPi, *BDTout, *is_signal), Import(*inTree_), Cut(DATA_selection));
// 	sidebands_data = (RooDataSet*)sidebands_data->reduce(*M_B0); //only keep *M_B0 var

// 	/*** SGN DATA ***/
	
// 	// retrieve MC data
// 	TString MCselection = SGNselection;
// 	MCselection.Append(Form(" && M_X3872 > %f && M_X3872 < %f", JpsiPiPi_SR_low, JpsiPiPi_SR_high)); //same cut on MC as on data

// 	RooDataSet* MC_data = new RooDataSet("MC_data", "MC_data", RooArgSet(*M_B0, *M_X3872, *M_K0s, *M_PiPi, *BDTout, *is_signal), Import(*inTree_), Cut(MCselection));
// 	MC_data = (RooDataSet*)MC_data->reduce(*M_B0); // reduce to only *M_B0 variable
	
// 	// // Rescale, i.e. set weight equal to McNormPerYear_X3872[year]
// 	// RooFormulaVar wFunc("w","event weight", Form("(%f)", McNormPerYear_X3872[year_]), *M_B0);
// 	// RooRealVar* w = (RooRealVar*)MC_data->addColumn(wFunc);

// 	// // redefine dataset with weightVar specified to rescale
// 	// MC_data = new RooDataSet(MC_data->GetName(), MC_data->GetTitle(), MC_data, *MC_data->get(), "", w->GetName());

// 	/*** TOTAL PSEUDODATASET ***/
// 	wsFit->import(*sidebands_data);
// 	wsFit->import(*MC_data);
// 	wsFit->import(*gen_data);
// 	wsFit->writeToFile("pseudo_ws.root");
	
// 	// RooDataSet* pseudodata = new RooDataSet("pseudodata", "pseudodata", RooArgSet(*M_B0), Import(*sidebands_data), Import(*MC_data), Import(*gen_data));
// 	// RooDataSet* pseudodata = (RooDataSet*)sidebands_data->Clone("pseudodata");
// 	// pseudodata->append(*MC_data);
// 	// pseudodata->append(*gen_data);

// 	RooDataHist* gen_data_binned = new RooDataHist("gen_data_binned", "gen_data_binned", *M_B0, *gen_data);
// 	RooDataHist* sidebands_data_binned = new RooDataHist("sidebands_data_binned", "sidebands_data_binned", *M_B0, *sidebands_data);

// 	RooDataHist* pseudodata = new RooDataHist("pseudodata", "pseudodata", *M_B0, *MC_data, McNormPerYear_X3872[year_]);
// 	pseudodata->add(*gen_data_binned, Form("(M_B0 > %f && M_B0 < %f)", B0_lSB_high, B0_rSB_low));
// 	pseudodata->add(*sidebands_data_binned, Form("(M_B0 < %f || M_B0 > %f)", B0_lSB_high, B0_rSB_low));

// 	// Debugging plots
// 	TCanvas* c = new TCanvas("c", "c", 600, 2400);
// 	c->Divide(1,4);
	
// 	c->cd(1);
// 	RooPlot* frame1 = M_B0->frame(Name("frame1"),Bins(55), Title("Sidebands data"));
// 	sidebands_data_binned->plotOn(frame1);
// 	frame1->Draw();

// 	c->cd(2);
// 	RooPlot* frame2 = M_B0->frame(Name("frame2"),Bins(55), Title("MC data"));
// 	MC_data->plotOn(frame2);
// 	frame2->Draw();

// 	c->cd(3);
// 	RooPlot* frame3 = M_B0->frame(Name("frame3"),Bins(55), Title("Generated data (between sidebands)"));
// 	gen_data_binned->plotOn(frame3);
// 	frame3->Draw();

// 	c->cd(4);
// 	RooPlot* frame4 = M_B0->frame(Name("frame4"),Bins(55), Title("Pseudodataset = sidebands + MC + generated"));
// 	pseudodata->plotOn(frame4);
// 	frame4->Draw();

// 	c->SaveAs(outPath_ + "pseudodataset.pdf");

// 	// save pseudodataset to workspace
// 	wsFit->import(*pseudodata);
// } // dataset


void MVAoptimizer::create_pseudodataset(){

	// retrieve all vars
	RooRealVar* M_B0 = (RooRealVar*) wsFit->var("M_B0");
	RooRealVar* M_X3872 = (RooRealVar*) wsFit->var("M_X3872");
	RooRealVar* M_K0s = (RooRealVar*) wsFit->var("M_K0s");
	RooRealVar* M_PiPi = (RooRealVar*) wsFit->var("M_PiPi");
	RooRealVar* BDTout = (RooRealVar*) wsFit->var("BDTout");
	RooRealVar* is_signal = new RooRealVar("is_signal", "Flag for MC or data", 0, 1, ""); //needed for selection

	// === RETRIEVE MODELS === //

	/*** BKG MODEL ***/

	// retrieve BKGmodel from workspace
	RooAbsPdf* BKGmodel = wsFit->pdf("BKGmodel"); //NOTE: after full model fit, BKGmodel params are updated to latest values

	// retrieve sidebands data (for plotting)
	TString sb_selection = DATAselection;
	sb_selection.Append(Form(" && (M_B0 < %f || M_B0 > %f)", B0_lSB_high, B0_rSB_low)); //select B0 sidebands
	sb_selection.Append(Form(" && M_X3872 > %f && M_X3872 < %f", JpsiPiPi_SR_low, JpsiPiPi_SR_high)); //select X (or Psi2S) SR

	std::cout << "checkpoint" << std::endl;
	RooDataSet* sidebands_data = new RooDataSet("sidebands_data", "sidebands_data", RooArgSet(*M_B0, *M_X3872, *M_K0s, *M_PiPi, *BDTout, *is_signal), Import(*inTree_), Cut(sb_selection));
	std::cout << "checkpoint" << std::endl;

	// retrieve n_bkg from full fit yield
	RooFitResult* ResBKGmodel = (RooFitResult*)wsFit->obj("fitresult_BKGmodel_data");
	double n_bkg_comb = ((RooRealVar*)(ResBKGmodel->floatParsFinal().find("n_comb")))->getVal();
	double n_bkg_jpsix = ((RooRealVar*)(ResBKGmodel->floatParsFinal().find("n_jpsix")))->getVal();
	double n_bkg = n_bkg_comb + n_bkg_jpsix;

	/*** SGN MODEL ***/
	
	RooAbsPdf* SGNmodel = wsFit->pdf("SGNmodel");

	// retrieve MC dataset to get number of events
	TString MCselection = SGNselection;
	MCselection.Append(Form(" && M_X3872 > %f && M_X3872 < %f", JpsiPiPi_SR_low, JpsiPiPi_SR_high)); //same cut on MC as on data
	RooDataSet* MC_data = new RooDataSet("MC_data", "MC_data", RooArgSet(*M_B0, *M_X3872, *M_K0s, *M_PiPi, *BDTout, *is_signal), Import(*inTree_), Cut(MCselection));
	MC_data = (RooDataSet*)MC_data->reduce(*M_B0); // reduce to only *M_B0 variable

	double n_sgn = McNormPerYear_X3872[year_] * MC_data->sumEntries();

	std::cout << "nsgn = " << n_sgn << " = " << McNormPerYear_X3872[year_] << " * " << MC_data->sumEntries() << std::endl;

	// === GENERATE DATA === //

	// bkg
	RooDataSet* gendata_bkg = BKGmodel->generate(*M_B0, NumEvents(n_bkg), Name("gendata_bkg"));
	RooDataHist* gendata_bkg_binned = new RooDataHist("gendata_bkg_binned", "gendata_bkg_binned", *M_B0, *gendata_bkg);

	// sgn
	RooDataSet* gendata_sgn = SGNmodel->generate(*M_B0, NumEvents((int)n_sgn), Name("gendata_sgn"));
	RooDataHist* gendata_sgn_binned = new RooDataHist("gendata_sgn_binned", "gendata_sgn_binned", *M_B0, *gendata_sgn);
	
	// sum
	RooDataHist* pseudodata = new RooDataHist("pseudodata", "pseudodata", *M_B0, *gendata_bkg_binned);
	pseudodata->add(*gendata_sgn_binned);


	// Debugging plots
	TCanvas* c = new TCanvas("c", "c", 600, 600*3);
	c->Divide(1,3);
	
	c->cd(1);
	RooPlot* frame1 = M_B0->frame(Name("frame1"),Bins(55), Title("Sidebands fit"));
	sidebands_data->plotOn(frame1, MarkerColor(kBlue), Name("sidebands_data_plot"));
	gendata_bkg->plotOn(frame1);
	BKGmodel->plotOn(frame1, Range("FULLregion"), NormRange("B0_lSB,B0_rSB"), LineColor(kRed), Name("bkgmodel_curve"));
	frame1->Draw();

	auto legend_sb = new TLegend();
	legend_sb->AddEntry("sidebands_data_plot", "Data", "ep");
	legend_sb->AddEntry("gendata_bkg", "Generated data", "ep");
	legend_sb->AddEntry("bkgmodel_curve", "Sidebands fit", "l");
	legend_sb->Draw();

	c->cd(2);
	RooPlot* frame2 = M_B0->frame(Name("frame2"),Bins(55), Title("MC fit"));
	gendata_sgn->plotOn(frame2);
	SGNmodel->plotOn(frame2, LineColor(kRed), Name("sgnmodel_curve"));
	frame2->Draw();

	auto legend_MC = new TLegend();
	legend_MC->AddEntry("gendata_sgn", "Generated MC", "ep");
	legend_MC->AddEntry("sgnmodel_curve", "MC fit", "l");
	legend_MC->Draw();

	c->cd(3);
	RooPlot* frame3 = M_B0->frame(Name("frame3"),Bins(55), Title("Pseudodataset"));
	pseudodata->plotOn(frame3);
	frame3->Draw();

	c->SaveAs(outPath_ + "pseudodataset.pdf");

	// Repeat plot separately for each box

	// Save plot data to .root file for external plots

	RooWorkspace *wout = new RooWorkspace("wout", "pseudodataset workspace");
	
	wout->import(*sidebands_data);
	wout->import(*gendata_bkg_binned);
	wout->import(*BKGmodel);
	wout->import(*MC_data);
	wout->import(*gendata_sgn_binned);
	wout->import(*SGNmodel);	
	wout->import(*pseudodata);

	wout->writeToFile("pseudodataset_plot_ws.root");
	
	// import data to workspace
	wsFit->import(*pseudodata);

}


// load S/B regions
void MVAoptimizer::set_SBregions(){

    // B0 sidebands limits
	// == from 4 to 10 sigma around mean [from MC fit]
    B0_lSB_low  = 5.1; //5.09783;
    B0_lSB_high = 5.18903;
    B0_rSB_low   = 5.37143;
    B0_rSB_high = 5.46262;

	if(channel_ == "X3872"){
		B0_lSB_low = 5;
		B0_rSB_high = 5.55;
	} 

    // B0 signal region
	// == 3 sigma around mean [from MC fit]
	B0_SR_low  = 5.21183;
    B0_SR_high = 5.34863;

    // JpsiPiPi signal region
    JpsiPiPi_SR_low = 3.82562;
    JpsiPiPi_SR_high= 3.91839;
    if(channel_ == "Psi2S"){
        JpsiPiPi_SR_low = 3.65831; 
        JpsiPiPi_SR_high= 3.71493; 
    }

    // K0s signal region
    K0s_SR_low  = 0.464637;
    K0s_SR_high = 0.530962;


    // plot limits
    MB_low = 5.0 , MB_high = 5.55;
    //MJpsiPiPi_low, MJpsiPiPi_high;
    //MK0s_low, MK0s_high;
    MPiPi_low = 0.3, MPiPi_high = 1.;
    BDTout_low = -10, BDTout_high = 5;


}//set_SBregions()

void MVAoptimizer::set_selection(){
    // selection on the K0s SR
    selection_ = Form("M_K0s > %f && M_K0s < %f && M_PiPi < 0.85 && M_X3872 < 4.0", K0s_SR_low, K0s_SR_high);
    // selection on the JpsiPiPi SR
    SRselection = selection_;
    // selection on the MVA 
    selection_.Append(Form(" && BDTout > %f && M_PiPi > %f", BDToutCut_, MpipiCut_)); 

    std::cout << " [SELECTION] to both data and MC " << selection_ << std::endl;
    SGNselection  = "is_signal && " + selection_;
    DATAselection = "!is_signal && " + selection_;

}// set_selection()


double MVAoptimizer::Nbkg_extraction(){

    double Nbkg;
    
	// ==== ROOFIT SET UP ==== //
	// prepare B0 mass dataset
	RooRealVar M_B0   ("M_B0", " B0 mass in data", MB_low, MB_high, "GeV");
    RooRealVar M_X3872("M_X3872", " JpsiPiPi mass in data", 3.4, 4.9, "GeV");
    RooRealVar M_K0s  ("M_K0s", " K0short mass in data", 0.4, 0.6, "GeV");
    RooRealVar M_PiPi ("M_PiPi", "PiPi mass", MPiPi_low, MPiPi_high, "GeV");
    RooRealVar BDTout ("BDTout", " BDTout in data", BDTout_low, BDTout_high, "");
	// define signal & sidebands region
	M_B0.setRange(  "B0_lSB"      , B0_lSB_low , B0_lSB_high);
	M_B0.setRange(  "B0_rSB"      , B0_rSB_low , B0_rSB_high);
	M_B0.setRange("BLINDregion"   , B0_lSB_high, B0_rSB_low);
	M_B0.setRange("FULLregion"    , MB_low,      MB_high);
	M_B0.setRange("SGNregion"     , B0_SR_low  , B0_SR_high);

    
	TString selection;
    RooDataSet tmp_data("data", "tmp_data", RooArgSet(M_B0, M_X3872, M_K0s, M_PiPi, BDTout), Import(*inTree_), Cut(selection_));
	tmp_data.Print();
	wsFit->import(tmp_data);

	selection = selection_;
	selection.Append(Form(" && (M_B0 < %f || M_B0 > %f)", B0_lSB_high, B0_rSB_low)); //select B0 sidebands
	selection.Append(Form(" && M_X3872 > %f && M_X3872 < %f", JpsiPiPi_SR_low, JpsiPiPi_SR_high)); //select X (or Psi2S) SR
    RooDataSet* data = (RooDataSet *)tmp_data.reduce(M_B0, selection);
    data->Print();
    
	// define BACKGROUND MODEL
    // Fermi 
	RooRealVar Slope("Slope", "", 10.0, 1.0 , 30.);
	RooRealVar Flex("Flex", "", 5., 3. , 5.25);
	RooRealVar C("C", "", 0.05, 0.0, 0.1);

    RooGenericPdf Fermi("Fermi", "", "1./(1. + exp(( @0 - @1)*@2)) + @3", RooArgList(M_B0, Flex, Slope, C));
    // Poisson
	RooRealVar SlopeP("SlopeP", "", -15.0, -35.0 , -5.);
	RooRealVar Th("Th", "", 4.75, 4.0 , 5.);
	RooRealVar Exp("Exp", "", 4.);
	
	RooGenericPdf Pois("Pois", "", "pow((@0 - @1), @4) * exp(( @0 - @1)*@2) + @3", RooArgList(M_B0, Th,  SlopeP, C, Exp));

    // Jpsi+X + combinatorial
    // exponential comb
    RooRealVar comb_coeff("comb_coeff","",-5,-15., -1);
    RooExponential pdf_comb("pdf_comb","",M_B0,comb_coeff);
    RooRealVar n_comb("n_comb","",80000,0.,1E6);
    // Jpsi+X (conjugated error function)
    RooRealVar jpsix_scale("jpsix_scale","",0.02,0.001,0.1);
    RooRealVar jpsix_shift("jpsix_shift","",5.15,5.12,5.16);
    RooGenericPdf pdf_jpsix("pdf_jpsix","","TMath::Erfc((@0-@1)/@2)",RooArgList(M_B0,jpsix_shift,jpsix_scale));
    RooRealVar n_jpsix("n_jpsix","",20000,0.,1E5);

    RooAddPdf BKGmodel("BKGmodel","",RooArgList(pdf_comb,pdf_jpsix),RooArgList(n_comb,n_jpsix));
    //RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", Fermi, bkg_yield);
	//RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", Pois, bkg_yield);


    // fit data
    //RooFitResult *ResBKGmodelSB2 = BKGmodel.fitTo(*data, Range("B0_rSB"), Save());
	//ResBKGmodelSB2->Print();
	//RooFitResult *ResBKGmodelSB1 = BKGmodel.fitTo(*data, Range("B0_lSB"), Save());
	//ResBKGmodelSB1->Print();
	RooFitResult *ResBKGmodel= BKGmodel.fitTo(*data, Range("B0_lSB,B0_rSB"), Save());
	ResBKGmodel->Print("v");
    wsFit->import(BKGmodel);
	wsFit->import(*ResBKGmodel);

    // plot results
    RooPlot* frame = M_B0.frame(Name("FitPlot"),Bins(55));
	data->plotOn(frame);
	BKGmodel.plotOn(frame, Range("FULLregion"), RooFit::NormRange("B0_lSB,B0_rSB"));
	BKGmodel.plotOn(frame,  Components(pdf_jpsix), LineStyle(kDashed), LineColor(kRed), Range("FULLregion"), RooFit::NormRange("B0_lSB,B0_rSB"));
	BKGmodel.plotOn(frame,  Components(pdf_comb), LineStyle(kDashed), LineColor(kGray), Range("FULLregion"), RooFit::NormRange("B0_lSB,B0_rSB"));
	BKGmodel.paramOn(frame, Layout(0.60));
	frame->SetTitle(Form("BDT out > %.2f Mpipi > %.3f", BDToutCut_, MpipiCut_));
	wsFit->import(*frame);

	// nFitParams = 2 normalizations + 1 exp coeff + 2 errf params
	TText *TxtChi2= new TText(5.2, frame->GetMaximum()*0.85, Form("Chi^2 = %.3f", frame->chiSquare(5)));
   	TxtChi2->SetTextSize(0.04);
   	TxtChi2->SetTextColor(kRed);
   	frame->addObject(TxtChi2);
	std::cout << " ---> Chi^2 = " << frame->chiSquare(5) << std::endl;

	TH2 *h_Corr = ResBKGmodel->correlationHist();
	
	// ==== Nbkg ==== //
	RooAbsReal* IntSreg = BKGmodel.createIntegral(M_B0, NormSet(M_B0), Range("SGNregion")); //integrate sgn region
	double Is = IntSreg->getVal();
	std::cout << " - Integral in signal region " << Is << std::endl;
	RooAbsReal* IntBreg = BKGmodel.createIntegral(M_B0, NormSet(M_B0), Range("B0_lSB,B0_rSB"));   //integrate bkg region
	double Ib = IntBreg->getVal();
	std::cout << " - Integral in sidebands region " << Ib << std::endl;

	//Nbkg = data->sumEntries() * Is/Ib; 
	Nbkg = inTree_->GetEntries(selection + Form("&& M_B0 > %f && M_B0 < %f", B0_lSB_low, B0_rSB_high))*Is/Ib;
	std::cout << " - Nbkg in signal region " << Nbkg << std::endl;

	// ==== SAVE ON FILE ==== //
	TString outFileName = "./outRoot/PunziOptimization/fit_results/FitResults_B0sidebands_" + year_+"_";
	if (BDToutCut_ < 0.) outFileName.Append(Form("m"));
	outFileName.Append(Form("%.0f_Mpipi%.0f", fabs(BDToutCut_)*100., MpipiCut_*1000) + channel_ + ".root");

	TFile* outFitFile = new TFile(outFileName, "RECREATE");
	frame->Write("FitPlot");
	ResBKGmodel->Write("FitResults");
	h_Corr->Write();
	wsFit->writeToFile("prova_workspace.root");

	outFitFile->Close();
	std::cout << " [OUT] written on " << outFileName << std::endl;
	//wsFit->Print();
    return Nbkg;

}//Nbkg_extraction

double MVAoptimizer::Nsgn_extraction(){

	// double Nsgn = inTree_->GetEntries(SGNselection);
	//std::cout << " - Nsgn in signal region " <<  Nsgn << std::endl;
	
	std::cout << "#### NSGN EXTRACTION " << std::endl;

	// ==== ROOFIT SET UP ==== //

	// retrieve fit var
	RooRealVar* M_B0 = wsFit->var("M_B0");
	RooRealVar* M_X3872 = wsFit->var("M_X3872");
	RooRealVar* M_K0s = wsFit->var("M_K0s");
	RooRealVar* M_PiPi = wsFit->var("M_PiPi");
	RooRealVar* BDTout = wsFit->var("BDTout");
	RooRealVar* is_signal = new RooRealVar("is_signal", "Flag for MC or data", 0, 1, ""); //needed for selection
	
	// retrieve data and only select MC events
	TString selection = SGNselection;
	selection.Append(Form(" && M_X3872 > %f && M_X3872 < %f", JpsiPiPi_SR_low, JpsiPiPi_SR_high));
	RooDataSet* tmpdata = new RooDataSet("tmpdata", "tmpdata", RooArgSet(*M_B0, *M_X3872, *M_K0s, *M_PiPi, *BDTout, *is_signal), Import(*inTree_), Cut(selection));
	RooDataSet* data = (RooDataSet*)tmpdata->reduce(RooArgSet(*M_B0)); // reduce to only M_B0 variable


	// Building SGN model

	// Crystal Ball
	// RooRealVar B0mass ("B0mass" , "", 5.279, 5.275, 5.285);
	// RooRealVar B0sigma_CB("B0sigma_CB", "", 0.007, 0.005, 0.1);
	// RooRealVar alpha("alpha", "", 1.5, 1.,  2.);
	// RooRealVar     N(    "N", "", 1., .1,  10.);
	RooRealVar B0mass ("B0mass" , "", 5.279, 5.275, 5.285);
	RooRealVar B0sigma_CB("B0sigma_CB", "", 0.007, 0.005, 0.1);
	RooRealVar alphaL("alphaL", "", 1., .5,  2.);
	RooRealVar alphaR("alphaR", "", 1., .5,  2.); //1.5 1 2
	RooRealVar     NR(    "NR", "", 5., 1.,  20.); //1. .1 10 [1 1 10]
	RooRealVar     NL(    "NL", "", 5., 1.,  50.); // 1. .1 10 [1 1 10]

	// RooCBShape CBsignal("CBsignal", "", *M_B0, B0mass, B0sigma_CB, alpha, N);
	// RooProdPdf SGNpdfCB("SGNpdfCB", "SGNpdfCB", CBsignal ); 
	RooCrystalBall DCBsignal("DCBsignal", "", *M_B0, B0mass, B0sigma_CB, alphaL, NL, alphaR, NR);

	// Gaussian
	RooRealVar B0sigma_G("B0sigma_G", "", 0.06, 0.001  , 0.1);

	RooGaussian Gsignal("Gsignal", "", *M_B0, B0mass, B0sigma_G);

	// CB + G
	// RooRealVar f("f", "", 0., 1.);
	// RooAddPdf SGNmodel("SGNmodel", "SGNmodel", RooArgList(CBsignal, Gsignal), f);
	RooProdPdf SGNmodel("SGNmodel", "SGNmodel", RooArgList(DCBsignal));

	// === FIT === //
	RooFitResult *ResSGNmodel= SGNmodel.fitTo(*data, Range("SGNregion"), Save());

	// save result to workspace
	ResSGNmodel->Print("v");
	wsFit->import(SGNmodel);
	wsFit->import(*ResSGNmodel);


	// ==== Nsgn ==== //

	RooAbsReal* IntSreg = SGNmodel.createIntegral(*M_B0, NormSet(*M_B0), Range("SGNregion")); //integrate sgn region
	double Nsgn = IntSreg->getVal();
	std::cout << " ### Nsgn in signal region " << Nsgn << std::endl;

	// === PLOT === //

	// Data + fit
	RooPlot* frame = M_B0->frame(Name("SGNFitPlot"),Bins(55));
	data->plotOn(frame, Name("mcpoints"));
	SGNmodel.plotOn(frame, Range("FULLregion"), RooFit::NormRange("SGNregion"), Name("sgnmodel_full"));
	// SGNmodel.plotOn(frame, Name("sgnmodel_cb"), Components(CBsignal), LineStyle(kDashed), LineColor(kRed), Range("FULLregion"), RooFit::NormRange("SGNregion"));
	// SGNmodel.plotOn(frame, Name("sgnmodel_g"), Components(Gsignal), LineStyle(kDashed), LineColor(kGray), Range("FULLregion"), RooFit::NormRange("SGNregion"));
	SGNmodel.paramOn(frame, Layout(0.6, 0.9, 0.9));
	frame->SetTitle(Form("MC SGN fit, BDT out > %.2f Mpipi > %.3f", BDToutCut_, MpipiCut_));
	wsFit->import(*frame);

	// chi2 
	TText *TxtChi2= new TText(5.05, frame->GetMaximum()*0.6, Form("Chi^2 = %.3f", frame->chiSquare("sgnmodel_full", "mcpoints", 6)));
	TxtChi2->SetTextSize(0.04);
	TxtChi2->SetTextColor(kRed);
	frame->addObject(TxtChi2);
	std::cout << " ---> Chi^2 = " << frame->chiSquare("sgnmodel_full", "mcpoints", 6) << std::endl;

	// canvas	
	TCanvas* c1 = new TCanvas();

	frame->Draw();

	// legend
	auto legend = new TLegend(0.15, 0.70,.45,.83);
	legend->SetTextSize(0.035);
	legend->SetBorderSize(0);
	legend->AddEntry("mcpoints", "MC", "lp");
	legend->AddEntry("sgnmodel_full", "Total fit", "l");
	// legend->AddEntry("sgnmodel_cb", "Crystal ball component", "l");
	// legend->AddEntry("sgnmodel_g", "Gaussian component", "l");
	legend->Draw();
	// gPad->BuildLegend();

	c1->SaveAs(outPath_ + "SGNfit_" + Form("X%.0f_Mpipi%.0f_%s_", fabs(BDToutCut_)*100., MpipiCut_*1000, year_.Data()) + channel_ + "_CBonly.pdf");
	
	Nsgn *= LumiNorm_;

	return Nsgn;

}//Nsgn_extraction()

double MVAoptimizer::EFFsgn_extraction(){

	double Ns = Nsgn_extraction();
	double Ntot = inTree_->GetEntries( "is_signal && " + SRselection)*LumiNorm_;

	return Ns/Ntot;

}//EFFsgn_extraction()

double MVAoptimizer::PunziSign(double* PSerr, double * Nbkg , double * Esgn ){

	double PunziS;
	const double b = 5.0;    // #sigmas corresp. to 5sigma significance level
	const double a = 2.0;    // #sigmas corresp. to CL (90%--> 1.2816) (95% --> 1.6448) (CMStwiki --> 2.)

	double B          = Nbkg_extraction();
    *Nbkg = B;
	double S          = Nsgn_extraction();
	double Seff       = EFFsgn_extraction();
    *Esgn = Seff;
	double Seff_error = 1./sqrt(S); // error square
	
	double Sign_denom        = b*b + 2.*a*sqrt(B) + b*sqrt(b*b + 4.*a*sqrt(B) + 4.*B );
	double Sign_denom_error  = a/sqrt(B) + b * (a / sqrt(B) + 2. ) / sqrt(b*b + 4.*a*sqrt(B) + 4.*B );	
	Sign_denom_error /= Sign_denom;	

	PunziS = Seff/Sign_denom;
	*PSerr = PunziS*sqrt(Seff_error + Sign_denom_error*Sign_denom_error);

	return PunziS;

}//PunziSign()

int MVAoptimizer::makeSGNvsBKGplot(){

	int Nbins = 55;
	double Mlow = MB_low, Mhigh = MB_high;
	TH1F* h_Data_B0 = new TH1F("Data_B0", "", Nbins, Mlow, Mhigh);
	TH1F* h_SGN_B0  = new TH1F("SGN_B0" , "", Nbins, Mlow, Mhigh); 
	double Mlow_X = 3.75, Mhigh_X = 4.0;
	if (channel_ == "Psi2S") { Mlow_X = 3.60, Mhigh_X = 3.75; }
	TH1F* h_Data_X = new TH1F("Data_X", "", 30, Mlow_X, Mhigh_X);
	TH1F* h_SGN_X  = new TH1F("SGN_X" , "", 30, Mlow_X, Mhigh_X);
	TH2F* h_Data_B0vsX = new TH2F("Data_B0vsX", "", Nbins, Mlow, Mhigh, 15, Mlow_X, Mhigh_X );
	TH2F* h_SGN_B0vsX  = new TH2F("SGN_B0vsX", "", Nbins, Mlow, Mhigh, 15, Mlow_X, Mhigh_X );

	double Mlow_K0s = K0s_SR_low*0.9, Mhigh_K0s = K0s_SR_high*1.1;
	TH1F* h_Data_K0s = new TH1F("Data_K0s", "", 30, Mlow_K0s, Mhigh_K0s);
	TH1F* h_SGN_K0s  = new TH1F("SGN_K0s" , "", 30, Mlow_K0s, Mhigh_K0s);

	TString BKG_selection = DATAselection;
	if (is_blind_) BKG_selection.Append(Form("&& (M_B0 < %f || M_B0 > %f)",  B0_lSB_high, B0_rSB_low));
	std::cout << "BKG_selection = " << BKG_selection << std::endl;
	TString SGN_selection = SGNselection;
    // data
	inTree_->Draw("M_B0>>Data_B0", BKG_selection + Form("&& M_X3872 > %f && M_X3872 < %f", JpsiPiPi_SR_low, JpsiPiPi_SR_high));
	inTree_->Draw("M_X3872>>Data_X", BKG_selection + Form("&& M_B0 > %f && M_B0 < %f", B0_SR_low, B0_SR_high));
	inTree_->Draw("M_K0s>>Data_K0s", BKG_selection + Form("&& M_B0 > %f && M_B0 < %f && M_X3872 > %f && M_X3872 < %f", B0_SR_low, B0_SR_high, JpsiPiPi_SR_low, JpsiPiPi_SR_high));
	inTree_->Draw("M_X3872:M_B0>>Data_B0vsX", BKG_selection);
    // MC
	inTree_->Draw("M_B0>>SGN_B0", SGN_selection + Form("&& M_X3872 > %f && M_X3872 < %f", JpsiPiPi_SR_low, JpsiPiPi_SR_high));
	std::cout << "\t### MC EVENTS IN SIGNAL REGION, M_B0 hist (pre-normalization): " << h_SGN_B0->Integral() << std::endl;
	if(channel_ == "Psi2S") h_SGN_B0->Scale(LumiNorm_);
	
	std::cout << "\t### MC EVENTS IN SIGNAL REGION, M_B0 hist (after normalization): " << h_SGN_B0->Integral() << std::endl;
	inTree_->Draw("M_X3872>>SGN_X", SGN_selection + Form("&& M_B0 > %f && M_B0 < %f", B0_SR_low, B0_SR_high));
	std::cout << "\t### MC EVENTS IN SIGNAL REGION, M_X3872 hist (pre-normalization): " << h_SGN_X->Integral() << std::endl;
	if(channel_ == "Psi2S") h_SGN_X->Scale(LumiNorm_);

	inTree_->Draw("M_K0s>>SGN_K0s", SGN_selection + Form("&& M_B0 > %f && M_B0 < %f && M_X3872 > %f && M_X3872 < %f", B0_SR_low, B0_SR_high, JpsiPiPi_SR_low, JpsiPiPi_SR_high));
	if(channel_ == "Psi2S") h_SGN_K0s->Scale(LumiNorm_);

	inTree_->Draw("M_X3872:M_B0>>SGN_B0vsX", SGN_selection);
	if(channel_ == "Psi2S") h_SGN_B0vsX->Scale(LumiNorm_);
    h_SGN_B0vsX->Add(h_Data_B0vsX);

    std::cout << " Jpsi PiPi under B0 peak >>> selection Data: " << BKG_selection << std::endl;
    std::cout << " Jpsi PiPi under B0 peak >>> DATA :" << h_Data_X->Integral() << std::endl;
    std::cout << " Jpsi PiPi under B0 peak >>> selection MC : " << SGN_selection << std::endl;
    std::cout << " Jpsi PiPi under B0 peak >>> MC :" << h_SGN_X->Integral() << std::endl;

	h_SGN_B0->GetXaxis()->SetTitle("M(B_{0}) [GeV]");
	h_SGN_B0->GetYaxis()->SetTitle(Form("Events/%.3f [GeV]", h_SGN_B0->GetXaxis()->GetBinWidth(1)));
	h_SGN_B0->GetXaxis()->SetTitleOffset(1.1); h_SGN_B0->GetXaxis()->SetTitleSize(0.04);
	h_SGN_B0->GetYaxis()->SetTitleSize(0.04);
	h_SGN_B0->SetLineWidth(3);
	h_SGN_B0->SetLineColor(kAzure + 1); h_SGN_B0->SetFillColorAlpha(kAzure + 1, 0.30);
	
	h_Data_B0->SetLineColor(kBlack);
	h_Data_B0->SetMarkerStyle(20);
	h_Data_B0->SetLineWidth(2);

	h_SGN_X->GetXaxis()->SetTitle("M(J#psi #pi #pi) [GeV]");
	h_SGN_X->GetYaxis()->SetTitle(Form("Events/%.3f [GeV]", h_SGN_X->GetXaxis()->GetBinWidth(1)));
	h_SGN_X->GetXaxis()->SetTitleOffset(1.1); h_SGN_X->GetXaxis()->SetTitleSize(0.04);
	h_SGN_X->GetYaxis()->SetTitleSize(0.04);
	h_SGN_X->SetLineWidth(3);
	h_SGN_X->SetLineColor(kViolet + 1); h_SGN_X->SetFillColorAlpha(kViolet + 1, 0.30);
	h_Data_X->SetLineColor(kBlack);
	h_Data_X->SetMarkerStyle(20);
	h_Data_X->SetLineWidth(2);

	h_SGN_K0s->GetXaxis()->SetTitle("M(K^{0}_{S}) [GeV]");
	h_SGN_K0s->GetYaxis()->SetTitle(Form("Events/%.3f [GeV]", h_SGN_K0s->GetXaxis()->GetBinWidth(1)));
	h_SGN_K0s->GetXaxis()->SetTitleOffset(1.1); h_SGN_K0s->GetXaxis()->SetTitleSize(0.04);
	h_SGN_K0s->GetYaxis()->SetTitleSize(0.04);
	h_SGN_K0s->SetLineWidth(3);
	h_SGN_K0s->SetLineColor(kGreen + 1); h_SGN_K0s->SetFillColorAlpha(kGreen + 1, 0.30);
	h_Data_K0s->SetLineColor(kBlack);
	h_Data_K0s->SetMarkerStyle(20);
	h_Data_K0s->SetLineWidth(2);

    h_SGN_B0vsX->GetXaxis()->SetTitle("M(B^{0}) [GeV]");
    h_SGN_B0vsX->GetYaxis()->SetTitle("M(J#psi #pi^{+}#pi^{-}) [GeV]");
    h_SGN_B0vsX->GetXaxis()->SetTitleOffset(1.1); h_SGN_B0vsX->GetXaxis()->SetTitleSize(0.04);
	h_SGN_B0vsX->GetYaxis()->SetTitleSize(0.04);
	h_SGN_B0vsX->SetLineWidth(3);
	h_SGN_B0vsX->SetLineColor(kAzure + 1); //h_SGN_B0vsX->SetFillColorAlpha(kAzure + 1, 0.30);
	h_Data_B0vsX->SetLineColor(kRed);
	h_Data_B0vsX->SetLineWidth(3);

	// get FIT
	RooCurve* FitCurve;
	if(channel_ == "Psi2S") FitCurve = (RooCurve*)((RooPlot*)wsFit->obj("FitPlot"))->getObject(1);//(RooCurve*)FitPlot->getObject(1);
	else FitCurve = (RooCurve*)((RooPlot*)wsFit->obj("FullFitPlot"))->findObject("fullmodelcurve");
	FitCurve->SetLineColor(kRed); FitCurve->SetLineStyle(0);

	// Legend 
	// auto legendB0 = new TLegend(0.53, 0.70,.83,.83);
	auto legendB0 = new TLegend(0.16, 0.70, 0.46, 0.83);
	legendB0->SetTextSize(0.035);
	legendB0->SetBorderSize(0);

	gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
	gStyle->SetOptStat(0);
	gStyle->SetLineWidth(3);

	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
	h_SGN_B0->GetYaxis()->SetRangeUser(0., 1.2 * std::max(h_SGN_B0->GetMaximum(), h_Data_B0->GetMaximum()));
	h_SGN_B0->Draw("HIST"); //HIST
	h_SGN_B0->Draw("E0 SAME"); //HIST
	legendB0->AddEntry(h_SGN_B0, "SIGNAL (MC)");
	legendB0->AddEntry(h_Data_B0,"DATA " + year_);
	// FitCurve->Draw("SAME");
	h_Data_B0->Draw("PE0 SAME");
	// legendB0->AddEntry(FitCurve, "SIDEBANDS-FIT");

	// insert line at SR borders
	TLine *l1 = new TLine(B0_SR_low, 0., B0_SR_low, 0.5 * std::max(h_SGN_B0->GetMaximum(), h_Data_B0->GetMaximum()));
	TLine *l2 = new TLine(B0_SR_high, 0., B0_SR_high, 0.5 * std::max(h_SGN_B0->GetMaximum(), h_Data_B0->GetMaximum()));
	l1->SetLineColor(kRed); l1->SetLineStyle(2); l1->SetLineWidth(2);
	l2->SetLineColor(kRed); l2->SetLineStyle(2); l2->SetLineWidth(2);	
	l1->Draw("SAME");
	l2->Draw("SAME");

	legendB0->Draw();
	//TString outPath = "./outRoot/PunziOptimization/";
	CMSxxx(c1);
	c1->SaveAs(outPath_ + "B0massCUT_" + Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000) + channel_ + ".png");
	c1->SaveAs(outPath_ + "B0massCUT_" + Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000.)+ channel_ + ".pdf");

	h_SGN_X->GetYaxis()->SetRangeUser(0., 1.2 * std::max(h_SGN_X->GetMaximum(), h_Data_X->GetMaximum()));
	h_SGN_X->Draw("HIST"); //HIST
	h_SGN_X->Draw("E0 SAME"); //HIST
	h_Data_X->Draw("PE0 SAME");
	auto legendX = new TLegend(0.16, 0.70, 0.46, 0.83);
	legendX->SetTextSize(0.035);
	legendX->SetBorderSize(0);

	legendX->AddEntry(h_SGN_X, "SIGNAL (MC)");
	legendX->AddEntry(h_Data_X, "DATA " + year_);
	legendX->Draw();
    //c1 = myRootLib::RatioPlot(h_Data_X,h_SGN_X);

	// insert line at SR borders
	TLine *l3 = new TLine(JpsiPiPi_SR_low, 0., JpsiPiPi_SR_low, 0.5 * std::max(h_SGN_X->GetMaximum(), h_Data_X->GetMaximum()));
	TLine *l4 = new TLine(JpsiPiPi_SR_high, 0., JpsiPiPi_SR_high, 0.5 * std::max(h_SGN_X->GetMaximum(), h_Data_X->GetMaximum()));
	l3->SetLineColor(kRed); l3->SetLineStyle(2); l3->SetLineWidth(2);
	l4->SetLineColor(kRed); l4->SetLineStyle(2); l4->SetLineWidth(2);
	l3->Draw("SAME");
	l4->Draw("SAME");

	CMSxxx(c1);
	c1->SaveAs(outPath_ + "JpsiPiPimassCUT_" + Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000) + channel_ + ".png");
	c1->SaveAs(outPath_ + "JpsiPiPimassCUT_" + Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000.)+ channel_ + ".pdf");

	h_SGN_K0s->GetYaxis()->SetRangeUser(0., 1.2 * std::max(h_SGN_K0s->GetMaximum(), h_Data_K0s->GetMaximum()));
	h_SGN_K0s->Draw("HIST"); //HIST
	h_SGN_K0s->Draw("E0 SAME"); //HIST
	h_Data_K0s->Draw("PE0 SAME");
	auto legendK0s = new TLegend(0.16, 0.70, 0.46, 0.83);
	legendK0s->SetTextSize(0.035);
	legendK0s->SetBorderSize(0);

	legendK0s->AddEntry(h_SGN_K0s, "SIGNAL (MC)");
	legendK0s->AddEntry(h_Data_K0s, "DATA " + year_);
	legendK0s->Draw();
    //c1 = myRootLib::RatioPlot(h_Data_K0s,h_SGN_K0s);

	// insert line at SR borders
	TLine *l5 = new TLine(K0s_SR_low, 0., K0s_SR_low, 0.5 * std::max(h_SGN_K0s->GetMaximum(), h_Data_K0s->GetMaximum()));
	TLine *l6 = new TLine(K0s_SR_high, 0., K0s_SR_high, 0.5 * std::max(h_SGN_K0s->GetMaximum(), h_Data_K0s->GetMaximum()));
	l5->SetLineColor(kRed); l5->SetLineStyle(2); l5->SetLineWidth(2);
	l6->SetLineColor(kRed); l6->SetLineStyle(2); l6->SetLineWidth(2);
	l5->Draw("SAME");
	l6->Draw("SAME");

	CMSxxx(c1);
	c1->SaveAs(outPath_ + "K0smassCUT_" + Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000) + channel_ + ".pdf");


    h_SGN_B0vsX->Draw("COLZ");
    //h_Data_B0vsX->Draw("COLZ same");
    gPad->RedrawAxis();
	CMSxxx(c1);
	c1->SaveAs(outPath_ + "B0vsX_massCUT_" + Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000)+ channel_ + ".png");
	c1->SaveAs(outPath_ + "B0vsX_massCUT_" + Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000)+ channel_ + ".pdf");
    

	//inFile->Close();
	delete h_Data_B0;
	delete h_SGN_B0;
	delete h_Data_X;
	delete h_SGN_X;
	delete h_Data_B0vsX;
	delete h_SGN_B0vsX;

	return 0;

}//makeSGNvsBKGplot()



void MVAoptimizer::Plot2D_BDToutMpipi(){

	double CorrF;
	const int BDTout_Nbins = 100; 
	//const double BDTout_low = -0.8, BDTout_high = 0.5;
	const int Mpipi_Nbins= 70; 
	//const double MPiPi_low = 0.2, MPiPi_high = .7;
	TH2F* h_CorrS = new TH2F("CorrS_BDTvsMpipi", "", Mpipi_Nbins, MPiPi_low, MPiPi_high, BDTout_Nbins, BDTout_low, BDTout_high);	
	TH2F* h_CorrD = new TH2F("CorrD_BDTvsMpipi", "", Mpipi_Nbins, MPiPi_low, MPiPi_high, BDTout_Nbins, BDTout_low, BDTout_high);	
	
	inTree_->Draw("BDTout:M_PiPi>>CorrS_BDTvsMpipi", "is_signal && " + SRselection);
	inTree_->Draw("BDTout:M_PiPi>>CorrD_BDTvsMpipi", "!is_signal && "+ SRselection);
	h_CorrS->GetXaxis()->SetTitle("M(#pi^{+} #pi^{-})(GeV)");
	h_CorrS->GetYaxis()->SetTitle("BDT output");
	h_CorrS->GetYaxis()->SetTitleOffset(1.06);
	h_CorrS->SetLineColor(kAzure +1);
	h_CorrS->SetFillColorAlpha(kAzure +1, 0.7);
	h_CorrD->SetLineColor(kRed);
	h_CorrD->SetFillColorAlpha(kRed, 0.50);
    TLine Mpipi_th  = TLine(MpipiCut_, BDTout_low, MpipiCut_, BDTout_high); 
    Mpipi_th.SetLineColor(13); Mpipi_th.SetLineStyle(7);
    TLine BDTout_th = TLine(MPiPi_low, BDToutCut_, MPiPi_high, BDToutCut_); 
    BDTout_th.SetLineColor(13); BDTout_th.SetLineStyle(7);

	// projection along BDTout and Mpipi
	// Mpipi
	TH1D* h_Mpipi_S = h_CorrS->ProjectionX();
	TH1D* h_Mpipi_B = h_CorrD->ProjectionX();
	h_Mpipi_B->Rebin(2.); h_Mpipi_S->Rebin(2.);
	h_Mpipi_B->Scale(1. / h_Mpipi_B->Integral()); h_Mpipi_S->Scale(1. / h_Mpipi_S->Integral());
	h_Mpipi_B->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV)");
	h_Mpipi_B->GetYaxis()->SetTitle(Form("Events / %.0f MeV", h_Mpipi_B->GetXaxis()->GetBinWidth(10)* 1000.));
	h_Mpipi_B->SetMaximum(1.2 * std::max(h_Mpipi_B->GetMaximum(), h_Mpipi_S->GetMaximum()));
	h_Mpipi_B->SetFillColor(0); h_Mpipi_S->SetFillColor(0);
	h_Mpipi_B->SetLineWidth(3); h_Mpipi_S->SetLineWidth(3);
	TH1D* h_BDTx_S = h_CorrS->ProjectionY();
	TH1D* h_BDTx_B = h_CorrD->ProjectionY();
	h_BDTx_B->Rebin(2.); h_BDTx_S->Rebin(2.);
	h_BDTx_B->Scale(1. / h_BDTx_B->Integral()); h_BDTx_S->Scale(1. / h_BDTx_S->Integral());
	h_BDTx_B->GetXaxis()->SetTitle("BDT output");
	h_BDTx_B->GetYaxis()->SetTitle(Form("Events / %.3f", h_BDTx_B->GetXaxis()->GetBinWidth(10)));
	h_BDTx_B->SetMaximum(1.2 * std::max(h_BDTx_B->GetMaximum(), h_BDTx_S->GetMaximum()));
	h_BDTx_B->SetFillColor(0); h_BDTx_S->SetFillColor(0);
	h_BDTx_B->SetLineWidth(3); h_BDTx_S->SetLineWidth(3);
	


	auto legend= new TLegend(0.525, 0.8,.89,.89);
	legend->SetTextSize(0.025);
	legend->SetBorderSize(0);
	
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);	
	gStyle->SetOptStat(0);
	legend->AddEntry(h_CorrS, Form("SIGNAL (MC) Corr = %.3f", h_CorrS->GetCorrelationFactor()));
	legend->AddEntry(h_CorrD, Form("DATA " + year_+ " Corr = %.3f", h_CorrD->GetCorrelationFactor()));
	h_CorrS->Draw("BOX ");
	h_CorrD->Draw("BOX SAME");
	h_CorrS->Draw("BOX SAME");
    Mpipi_th.Draw("same"); BDTout_th.Draw("same");
	legend->Draw();
    CMSxxx(c1);

	c1->SaveAs(outPath_ + "BDToutVSMpipi_UL" +  year_+ "_" + channel_+ ".png"); 
	c1->SaveAs(outPath_ + "BDToutVSMpipi_UL" +  year_+ "_" + channel_+ ".pdf"); 

	h_Mpipi_B->Draw("HIST"); h_Mpipi_S->Draw("HIST SAME");
	gPad->SetLeftMargin(0.15);
	legend->Draw();
	c1->SaveAs(outPath_ + "Mpipi_dataVSmc_UL" +  year_+ "_" + channel_+ ".png"); 
	c1->SaveAs(outPath_ + "Mpipi_dataVSmc_UL" +  year_+ "_" + channel_+ ".pdf"); 

	h_BDTx_B->Draw("HIST"); h_BDTx_S->Draw("HIST SAME");
	gPad->SetLeftMargin(0.15);
	legend->Draw();
	c1->SaveAs(outPath_ + "BDTout_dataVSmc_UL" +  year_+ "_" + channel_+ ".png"); 
	c1->SaveAs(outPath_ + "BDTout_dataVSmc_UL" +  year_+ "_" + channel_+ ".pdf"); 

}

RooRealVar MVAoptimizer::TotalFit_Psi2S(){

	// ==== ROOFIT SET UP ==== //
	// define selection
	TString DATA_selection = DATAselection;
	DATA_selection.Append(Form("&& M_X3872 > %f && M_X3872 < %f", JpsiPiPi_SR_low, JpsiPiPi_SR_high));
	if(is_blind_) DATA_selection.Append(Form(" && (M_B0 < %f || M_B0 > %f)", B0_lSB_high, B0_rSB_low)); //only select sidebands 

	RooRealVar* M_B0 = wsFit->var("M_B0");
	TH1F hh("hh", "", 55, MB_low, MB_high);
    inTree_->Draw("M_B0>>hh", DATA_selection);
    RooDataHist h_data("h_data", "h_data", *M_B0, Import(hh)); // <--- [NOTE!] Binned fit

    // take BKG model and  freeze Jpsi+X bkg
    RooAbsPdf* BKGmodel = wsFit->pdf("BKGmodel"); 
    // RooRealVar* jpsix_scale = wsFit->var("jpsix_scale"); jpsix_scale->setConstant(true);
    // RooRealVar* jpsix_shift = wsFit->var("jpsix_shift"); jpsix_shift->setConstant(true);
    // RooRealVar* n_jpsix = wsFit->var("n_jpsix"); n_jpsix->setConstant(true);
	RooRealVar* jpsix_scale = wsFit->var("jpsix_scale");
	jpsix_scale->setMin(jpsix_scale->getVal() - jpsix_scale->getError()*3); jpsix_scale->setMax(jpsix_scale->getVal() + jpsix_scale->getError()*3);
    RooRealVar* jpsix_shift = wsFit->var("jpsix_shift");
	jpsix_shift->setMin(jpsix_shift->getVal() - jpsix_shift->getError()*3); jpsix_shift->setMax(jpsix_shift->getVal() + jpsix_shift->getError()*3);
    RooRealVar* n_jpsix = wsFit->var("n_jpsix");
	n_jpsix->setMin(n_jpsix->getVal() - n_jpsix->getError()*3); n_jpsix->setMax(n_jpsix->getVal() + n_jpsix->getError()*3);

	RooRealVar* n_comb = wsFit->var("n_comb");
	n_comb->setMin(n_comb->getVal() - n_comb->getError()*3);
	n_comb->setMax(n_comb->getVal() + n_comb->getError()*3);
	RooRealVar* comb_coeff = wsFit->var("comb_coeff");
	comb_coeff->setMin(comb_coeff->getVal() - comb_coeff->getError()*3);
	comb_coeff->setMax(comb_coeff->getVal() + comb_coeff->getError()*3);

    // fit variables 
    // SIGNAL model CB + G 
    // Crystal Ball

	RooAbsPdf* SGNmodel = wsFit->pdf("SGNmodel");

	RooRealVar* B0mass = wsFit->var("B0mass");
	B0mass->setMin(5.279); B0mass->setMax(5.281);
	// B0mass->setMin(B0mass->getVal() - B0mass->getError());
	// B0mass->setMax(B0mass->getVal() + B0mass->getError());
	// std::cout << Form("### TOTAL FIT B0_mass range = [%f, %f]", B0mass->getMin(), B0mass->getMax()) << std::endl;
	RooRealVar* B0sigma = wsFit->var("B0sigma_CB"); 

	//initialize other vars to MC fit value + restrict range
	RooRealVar* alphaL = wsFit->var("alphaL");
	alphaL->setMin(alphaL->getVal() - alphaL->getError()*3);
	alphaL->setMax(alphaL->getVal() + alphaL->getError()*3);
	RooRealVar* alphaR = wsFit->var("alphaR");
	alphaR->setMin(alphaR->getVal() - alphaR->getError()*3);
	alphaR->setMax(alphaR->getVal() + alphaR->getError()*3);
	RooRealVar* NL = wsFit->var("NL");
	NL->setMin(NL->getVal() - NL->getError()*3);
	NL->setMax(NL->getVal() + NL->getError()*3);
	RooRealVar* NR = wsFit->var("NR");
	NR->setMin(NR->getVal() - NR->getError()*3);
	NR->setMax(NR->getVal() + NR->getError()*3);

    // FULL model
    RooRealVar sgn_yield("sgn_yield", "", 100, 50, 1E6);
    RooRealVar bkg_yield("bkg_yield", "", 100, 80, 1E6);
    RooAddPdf FULLmodel("FULLmodel", "FULLmodel", RooArgList(*BKGmodel, *SGNmodel), RooArgList(bkg_yield,sgn_yield));

    // FIT
    RooFitResult *ResFULLmodel = FULLmodel.fitTo(h_data, Extended(kTRUE), Range("FULLregion"), Save());
    ResFULLmodel->Print("v");
    
    // DRAW RESULTS
    RooPlot* frame = M_B0->frame(Name("FullFitPlot"), Bins(55), Title(""));
    h_data.plotOn(frame, Name("h_data"));
    FULLmodel.plotOn(frame, Components(*SGNmodel), LineColor(kGreen), LineStyle(kDashed), NormRange("FULLregion"), Name("fullmodel_sgn"));
    FULLmodel.plotOn(frame, Components(*BKGmodel), LineColor(kRed), LineWidth(2), LineStyle(kDashed), NormRange("FULLregion"), Name("fullmodel_bkg"));
    FULLmodel.plotOn(frame, Name("fullmodelcurve"), NormRange("FULLregion"));
	// FULLmodel.paramOn(frame, Layout(0.55, 0.9, 0.95));
    RooHist *hpull = frame->pullHist("h_data", "fullmodelcurve", false);

    // chi-square
    TText *TxtChi2= new TText(5.1, frame->GetMaximum()*0.85, Form("ChiSq = %.3f", frame->chiSquare(13)));
    TxtChi2->SetTextSize(0.035);   TxtChi2->SetTextColor(kRed);
    // frame->addObject(TxtChi2);
    std::cout << " ----> Chi2 SIGNAL FIT \t" << frame->chiSquare(13) << std::endl;

    // plot pull
    RooPlot* frame_pull = M_B0->frame(Bins(55), Title("Pulls"));
    frame_pull->addPlotable(hpull, "P");

	// Save to WS
	wsFit->import(*frame);
	wsFit->import(FULLmodel);
	wsFit->import(*ResFULLmodel);
	wsFit->import(*frame_pull);

    //SGN and BKG #events
    double I_SGN_T  = SGNmodel->createIntegral(*M_B0, NormSet(*M_B0), Range("FULLregion"))->getVal(); //integrate sgn region
    double I_SGN_Sr = SGNmodel->createIntegral(*M_B0, NormSet(*M_B0), Range("SGNregion"))->getVal(); //integrate sgn region
    double nSGN = ((RooRealVar*)(ResFULLmodel->floatParsFinal().find("sgn_yield")))->getVal();
	double nSGN_error = ((RooRealVar*)(ResFULLmodel->floatParsFinal().find("sgn_yield")))->getError();
    double numSignal_SR = nSGN*I_SGN_Sr/I_SGN_T;
	double numSignal_SR_error = nSGN_error*I_SGN_Sr/I_SGN_T; //ignoring integral error propagation
    std::cout << Form(" [RESULT] nSGN = %.2f +- %.2f \t sgn-integral (SR)/(FULL) = %.2f / %.2f  --->> Nsignal = %.2f ", nSGN, nSGN_error, I_SGN_Sr, I_SGN_T, numSignal_SR) << std::endl;
    double I_BKG_T  = BKGmodel->createIntegral(*M_B0, NormSet(*M_B0), Range("FULLregion"))->getVal(); //integrate sgn region
    double I_BKG_Sr = BKGmodel->createIntegral(*M_B0, NormSet(*M_B0), Range("SGNregion"))->getVal(); //integrate sgn region
    double nBKG = ((RooRealVar*)(ResFULLmodel->floatParsFinal().find("bkg_yield")))->getVal();
    double numBackground_SR = nBKG*I_BKG_Sr/I_BKG_T;
    std::cout << Form(" [RESULT] nBKG = %.2f \t bkg-integral (SR)/(FULL) = %.2f / %.2f  --->> Nbackground = %.2f ", nBKG, I_BKG_Sr, I_BKG_T, numBackground_SR) << std::endl;
    TText *TxtNs = new TText(5.1, frame->GetMaximum()*0.8, Form("N SGN = %.2f ", numSignal_SR ));
    TText *TxtNb = new TText(5.1, frame->GetMaximum()*0.75, Form("N BKG = %.2f ", numBackground_SR));
    // frame->addObject(TxtNs);
    // frame->addObject(TxtNb);

    // TCanvas* cf = new TCanvas("cf","data fit", 500, 1000);
	TCanvas* cf = new TCanvas("cf","data fit", 500, 500);
    // cf->Divide(1,2);
    cf->cd(1);
    frame->Draw();

	gStyle->SetOptStat(0);

	// TH1F* h_SGN_B0  = new TH1F("SGN_B0" , "", 55, MB_low, MB_high); 
	// inTree_->Draw("M_B0>>SGN_B0", SGNselection + Form("&& M_X3872 > %f && M_X3872 < %f", JpsiPiPi_SR_low, JpsiPiPi_SR_high), "goff");
	// h_SGN_B0->Scale(LumiNorm_);
	// h_SGN_B0->SetLineWidth(3);
	// h_SGN_B0->SetLineColor(33); h_SGN_B0->SetFillColorAlpha(33, 0.30);
	// h_SGN_B0->Draw("HIST SAME"); //HIST
	// h_SGN_B0->Draw("E0 SAME"); //HIST
	// frame->addObject(h_SGN_B0, "HIST SAME");
	// frame->addObject(h_SGN_B0, "E0 SAME");

	frame->GetXaxis()->SetTitle("m_{B^{0}} [GeV]");
	frame->Draw("SAME");

	auto legend = new TLegend(0.2, 0.70, .5, .83);
	legend->SetTextSize(0.035);
	legend->SetBorderSize(0);
	legend->AddEntry("h_data", "Data 2022", "lp");
	// legend->AddEntry("SGN_B0", "MC", "lp");
	legend->AddEntry("fullmodelcurve", "Total fit", "l");
	legend->AddEntry("fullmodel_sgn", "Signal comp.", "l");
	legend->AddEntry("fullmodel_bkg", "Background comp.", "l");
	legend->Draw();

	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.05);

    // cf->cd(2);
    // frame_pull->Draw();

	CMSxxx(cf);

    cf->SaveAs(outPath_ + "FullFit_" +Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000)+ year_ + "_" + channel_ +".png");
    cf->SaveAs(outPath_ + "FullFit_" + Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000)+ year_ + "_" + channel_ +".pdf");

	// save log-scale plot
	TCanvas* cf_log = new TCanvas("cf_log", "data fit, log scale", 800, 800);

	RooPlot* frame_log = M_B0->frame(Name("FullFitPlot"), Bins(55), Title(""));
    h_data.plotOn(frame_log, Name("h_data"));
    FULLmodel.plotOn(frame_log, Name("fullmodelcurve"), NormRange("FULLregion"));

	frame_log->Draw();
	gPad->SetLogy();
	// frame_log->SetMinimum(100);

	cf_log->SaveAs(outPath_ + "FullFit_" +Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000)+ year_ + "_" + channel_ +"_logY.pdf");

	// save ws to file
	wsFit->writeToFile("full_fit_ws.root");

	// return nsgn with error
	RooRealVar res_nsgn("res_nsgn", "Number of signal events in SR region, from fit", numSignal_SR);
	res_nsgn.setError(numSignal_SR_error);

    return res_nsgn;
	
}//TotalFit_Psi2S()	

RooRealVar MVAoptimizer::TotalFit_X3872(){

	// TotalFit_Psi2S();		// run total fit
	create_pseudodataset(); // imports "pseudodataset" RooDataSet to workspace

	// retrieve pseudodataset from workspace
	// RooDataSet* pseudodata = (RooDataSet*)wsFit->data("pseudodata");
	RooDataHist* pseudodata = (RooDataHist*)wsFit->data("pseudodata");

	// build FULLmodel [same as TotalFit_Psi2S]

    // // take BKG model and  freeze Jpsi+X bkg
    RooAbsPdf* BKGmodel = wsFit->pdf("BKGmodel"); 
	RooRealVar* jpsix_scale = wsFit->var("jpsix_scale");
	jpsix_scale->setMin(jpsix_scale->getVal() - jpsix_scale->getError()*1); jpsix_scale->setMax(jpsix_scale->getVal() + jpsix_scale->getError()*1);
    RooRealVar* jpsix_shift = wsFit->var("jpsix_shift");
	jpsix_shift->setMin(jpsix_shift->getVal() - jpsix_shift->getError()*1); jpsix_shift->setMax(jpsix_shift->getVal() + jpsix_shift->getError()*1);
    RooRealVar* n_jpsix = wsFit->var("n_jpsix");
	n_jpsix->setMin(n_jpsix->getVal() - n_jpsix->getError()*1); n_jpsix->setMax(n_jpsix->getVal() + n_jpsix->getError()*1);

	RooRealVar* n_comb = wsFit->var("n_comb");
	n_comb->setMin(n_comb->getVal() - n_comb->getError()*1);
	n_comb->setMax(n_comb->getVal() + n_comb->getError()*1);
	RooRealVar* comb_coeff = wsFit->var("comb_coeff");
	comb_coeff->setMin(comb_coeff->getVal() - comb_coeff->getError()*1);
	comb_coeff->setMax(comb_coeff->getVal() + comb_coeff->getError()*1);

    // fit variables 
    // SIGNAL model CB + G 
    // Crystal Ball

	RooAbsPdf* SGNmodel = wsFit->pdf("SGNmodel");

	RooRealVar* B0mass = wsFit->var("B0mass");
	// B0mass->setConstant();
	// B0mass->setMin(5.279); B0mass->setMax(5.281);
	B0mass->setMin(B0mass->getVal() - B0mass->getError()*1);
	B0mass->setMax(B0mass->getVal() + B0mass->getError()*1);
	// // std::cout << Form("### TOTAL FIT B0_mass range = [%f, %f]", B0mass->getMin(), B0mass->getMax()) << std::endl;
	RooRealVar* B0sigma = wsFit->var("B0sigma_CB"); 
	// B0sigma->setConstant();
	B0sigma->setMin(B0sigma->getVal() - B0sigma->getError()*1);
	B0sigma->setMax(B0sigma->getVal() + B0sigma->getError()*1);

	//initialize other vars to MC fit value + restrict range
	RooRealVar* alphaL = wsFit->var("alphaL");
	// alphaL->setConstant();
	alphaL->setMin(alphaL->getVal() - alphaL->getError()*1);
	alphaL->setMax(alphaL->getVal() + alphaL->getError()*1);
	RooRealVar* alphaR = wsFit->var("alphaR");
	// alphaR->setConstant();
	alphaR->setMin(alphaR->getVal() - alphaR->getError()*1);
	alphaR->setMax(alphaR->getVal() + alphaR->getError()*1);
	RooRealVar* NL = wsFit->var("NL");
	// NL->setConstant();
	NL->setMin(NL->getVal() - NL->getError()*1);
	NL->setMax(NL->getVal() + NL->getError()*1);
	RooRealVar* NR = wsFit->var("NR");
	// NR->setConstant();
	NR->setMin(NR->getVal() - NR->getError()*1);
	NR->setMax(NR->getVal() + NR->getError()*1);

    // FULL model
    RooRealVar sgn_yield("sgn_yield", "", 100, 50, 1E6);
    RooRealVar bkg_yield("bkg_yield", "", 100, 80, 1E6);
    RooAddPdf* FULLmodel = new RooAddPdf("FULLmodel", "FULLmodel", RooArgList(*BKGmodel, *SGNmodel), RooArgList(bkg_yield,sgn_yield));

	// fit
	RooFitResult *ResFULLmodel_pseudo = FULLmodel->fitTo(*pseudodata, Extended(kTRUE), Range("FULLregion"), Save());
    ResFULLmodel_pseudo->Print("v");

	// extract n_sgn in SR region
	RooRealVar* M_B0 = wsFit->var("M_B0");

	double I_SGN_T  = SGNmodel->createIntegral(*M_B0, NormSet(*M_B0), Range("FULLregion"))->getVal(); //integrate sgn region
    double I_SGN_Sr = SGNmodel->createIntegral(*M_B0, NormSet(*M_B0), Range("SGNregion"))->getVal(); //integrate sgn region
    double nSGN = ((RooRealVar*)(ResFULLmodel_pseudo->floatParsFinal().find("sgn_yield")))->getVal();
	double nSGN_error = ((RooRealVar*)(ResFULLmodel_pseudo->floatParsFinal().find("sgn_yield")))->getError();
    double numSignal_SR = nSGN*I_SGN_Sr/I_SGN_T;
	double numSignal_SR_error = nSGN_error*I_SGN_Sr/I_SGN_T; //ignoring integral error propagation

	// Draw results
	RooPlot* frame = M_B0->frame(Name("FullFitPlot"), Bins(55), Title(" "));
	frame->GetXaxis()->SetTitle("m_{B^{0}} [GeV]");
	pseudodata->plotOn(frame, Name("h_data"));
	FULLmodel->plotOn(frame, Components("SGNmodel"), LineColor(kGreen), LineStyle(kDashed), NormRange("FULLregion"), Name("fullmodel_sgn"));
	FULLmodel->plotOn(frame, Components("BKGmodel"), LineColor(kRed), LineWidth(2), LineStyle(kDashed), NormRange("FULLregion"), Name("fullmodel_bkg"));
	FULLmodel->plotOn(frame, Name("fullmodelcurve"), NormRange("FULLregion"));
	// FULLmodel->paramOn(frame, Layout(0.55, 0.9, 0.95));
	// frame->getAttText()->SetTextSize(0.02);

	//chi2
	TText *TxtChi2= new TText(5.1, frame->GetMaximum()*0.85, Form("Chi^2 = %.3f", frame->chiSquare(13)));
   	TxtChi2->SetTextSize(0.04);
   	TxtChi2->SetTextColor(kRed);
   	// frame->addObject(TxtChi2);
	std::cout << " ---> Chi^2 = " << frame->chiSquare(13) << std::endl;

	//n_sgn
	TText *TxtNs = new TText(5.1, frame->GetMaximum()*0.8, Form("N SGN = %.2f ", numSignal_SR));
    // frame->addObject(TxtNs);

	TCanvas* cpseudo = new TCanvas("cpseudo","pseudodata fit", 600, 600);
	CMSxxx(cpseudo, "simulation");
	frame->Draw();

	//legend
	auto legend = new TLegend(0.5, 0.70, .89, .89);
	legend->SetTextSize(0.035);
	legend->SetBorderSize(0);
	legend->AddEntry("h_data", "Pseudo-data", "lp");
	legend->AddEntry("fullmodelcurve", "Total fit", "l");
	legend->AddEntry("fullmodel_sgn", "Signal comp.", "l");
	legend->AddEntry("fullmodel_bkg", "Background comp.", "l");
	legend->Draw();

	CMSxxx(cpseudo);

	cpseudo->SaveAs(outPath_ + "FullFit_pseudodata_" +Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000)+ year_ + "_" + channel_ + ".pdf");

	// export result
	wsFit->import(*frame);
	wsFit->import(*FULLmodel);
	wsFit->import(*ResFULLmodel_pseudo);

	// return N_SGN with error
	RooRealVar res_nsgn_pseudo("res_nsgn_pseudo", "Number of signal events in SR region in pseudodataset, from fit", numSignal_SR);
	res_nsgn_pseudo.setError(numSignal_SR_error);
	return res_nsgn_pseudo;
}

void MVAoptimizer::CMSxxx(TCanvas* c, TString data_or_sim){
	c->cd();
	TLatex RunDetails; RunDetails.SetNDC(); 
	RunDetails.SetTextFont(61);
	RunDetails.SetTextAlign(10);
	RunDetails.SetTextSize(0.035);
	RunDetails.DrawLatex(.14, .92, "CMS " + data_or_sim);
	RunDetails.SetTextFont(52);
	RunDetails.SetTextSize(0.035);

	double x_min = 0.48;
	if(!data_or_sim.CompareTo("data")) x_min = 0.3;
	else if(!data_or_sim.CompareTo("simulation")) x_min = 0.39;

	RunDetails.DrawLatex(x_min, .92, "Private work");
	RunDetails.SetTextFont(42);
	RunDetails.SetTextSize(0.030);
	RunDetails.DrawLatex(.70, .91, Form("%.2f fb^{-1} (13.6 TeV)",LumiPerYear_Run3[year_]));

	gPad->SetLeftMargin(0.14);

}
