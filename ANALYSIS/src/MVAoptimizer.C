#include "../include/MVAoptimizer.h"

MVAoptimizer::MVAoptimizer(const TString& input, const float& BDTcut,const float& Mcut,const int& year, const TString & channel){

    year_ = year;
	LumiNorm_ = 0.025;
    channel_ = channel;
	std::cout<< Form(" Analyze %d dataset for " + channel_ + "channel -> LumiNorm %.4f", year_, LumiNorm_)<<std::endl;

    BDToutCut_ = BDTcut;
    MpipiCut_  = Mcut;
	is_blind_ = true;

    inPath_ = input;
    if(!gSystem->AccessPathName(inPath_)){
        inFile_ = TFile::Open(inPath_);
        inTree_ = (TTree*)inFile_->Get("CVtrainig_UL_2017");
        std::cout << " [+] loaded TTree from " << inPath_ << " with "<< inTree_->GetEntries() << " entries" << std::endl;
    }else std::cout << " ERROR : cannot open input file " << inPath_ << std::endl;

    set_SBregions();
    set_selection();

    outPath_ = Form("/eos/user/c/cbasile/www/B0toX3872K0s/MVA/CVtraining_UL_%d/", year_);

}// MVAoptimizer()

MVAoptimizer::~MVAoptimizer(){
    inFile_->Close();
	wsFit->Delete();
}//~MVAoptimizer()

// load S/B regions
void MVAoptimizer::set_SBregions(){

    // B0 sidebands limits
    B0_lSB_low  = 5.09783;
    B0_lSB_high = 5.18903;
    B0_rSB_low  = 5.37143;
    B0_rSB_high = 5.46262;

    // B0 signal region
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
    MB_low = 5.05 , MB_high = 5.50;
    //MJpsiPiPi_low, MJpsiPiPi_high;
    //MK0s_low, MK0s_high;
    MPiPi_low = 0.3, MPiPi_high = 1.;
    BDTout_low = -10, BDTout_high = 5;


}//set_SBregions()

void MVAoptimizer::set_selection(){
    // selection on the K0s SR
    selection_ = Form("M_K0s > %f && M_K0s < %f", K0s_SR_low, K0s_SR_high);
    // selection on the JpsiPiPi SR
    selection_.Append(Form(" && M_X3872 > %f && M_X3872 < %f", JpsiPiPi_SR_low, JpsiPiPi_SR_high));
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
	//double Mlow = 5.0, Mhigh = 5.6;	
	RooRealVar M_B0   ("M_B0", " B0 mass in data", MB_low, MB_high, "GeV");
    RooRealVar M_X3872("M_X3872", " JpsiPiPi mass in data", 3.4, 4.9, "GeV");
    RooRealVar M_K0s  ("M_K0s", " K0short mass in data", 0.4, 0.6, "GeV");
    RooRealVar M_PiPi ("M_PiPi", "PiPi mass", MPiPi_low, MPiPi_high, "GeV");
    RooRealVar BDTout ("BDTout", " BDTout in data", BDTout_low, BDTout_high, "");
	// define signal & sidebands region
	M_B0.setRange(  "B0_lSB"      , MB_low , B0_lSB_high);
	M_B0.setRange(  "B0_rSB"      , B0_rSB_low , MB_high);
	M_B0.setRange("BLINDregion"   , B0_lSB_high, B0_rSB_low);
	M_B0.setRange("FULLregion"    , MB_low,      MB_high);
	M_B0.setRange("SGNregion"     , B0_SR_low  , B0_SR_high);
	// cuts on the intermediate resonances
	TString selection;
	//if (channel_ == "X3872") selection = "M_X3872 > 3.75 && M_PiPi < 0.9";
	//else if (channel_ == "Psi2S") selection = "M_X3872 < 3.75";
    RooDataSet tmp_data("data", "tmp_data", RooArgSet(M_B0, M_X3872, M_K0s, M_PiPi, BDTout), Import(*inTree_), Cut(selection));
	tmp_data.Print();

	selection = selection_;
	selection.Append(Form(" && (M_B0 < %f || M_B0 > %f)", B0_lSB_high, B0_rSB_low));
    RooDataSet* data = (RooDataSet *)tmp_data.reduce(M_B0, selection);
    data->Print();
	wsFit->import(*data);
    
	// define BACKGROUND MODEL
    // Fermi 
	RooRealVar Slope("Slope", "", 10.0, 1.0 , 200.);
	RooRealVar Flex("Flex", "", 4., 1. , 7.);
	RooRealVar C("C", "", 0.5, 0. , 1.);

    RooGenericPdf Fermi("Fermi", "", "1./(1. + exp(( @0 - @1)*@2)) + @3", RooArgList(M_B0, Flex, Slope, C));
    // Poisson
	RooRealVar SlopeP("SlopeP", "", -15.0, -35.0 , -5.);
	RooRealVar Th("Th", "", 4.75, 4.5 , 5.);
	RooRealVar Exp("Exp", "", 4.);
	
	RooGenericPdf Pois("Pois", "", "pow((@0 - @1), @4) * exp(( @0 - @1)*@2) + @3", RooArgList(M_B0, Th,  SlopeP, C, Exp));
	
	
	RooRealVar bkg_yield("nBKG", "", 3000, 100, 100000);
    //RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", Fermi, bkg_yield);
	RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", Pois, bkg_yield);
    wsFit->import(BKGmodel);


    // fit data
    //RooFitResult *ResBKGmodelSB2 = BKGmodel.fitTo(*data, Range("B0_rSB"), Save());
	//ResBKGmodelSB2->Print();
	//RooFitResult *ResBKGmodelSB1 = BKGmodel.fitTo(*data, Range("B0_lSB"), Save());
	//ResBKGmodelSB1->Print();
	RooFitResult *ResBKGmodel= BKGmodel.fitTo(*data, Range("B0_lSB,B0_rSB"), Save());
	ResBKGmodel->Print("v");
	wsFit->import(*ResBKGmodel);

    // plot results
    RooPlot* frame = M_B0.frame(Name("FitPlot"),Bins(45));
	data->plotOn(frame);
	BKGmodel.plotOn(frame, Range("FULLregion"), RooFit::NormRange("B0_lSB,B0_rSB"));
	BKGmodel.paramOn(frame, Layout(0.60));
	frame->SetTitle(Form("BDT out > %.2f Mpipi > %.3f", BDToutCut_, MpipiCut_));
	wsFit->import(*frame);

	
	TText *TxtChi2= new TText(5.2, frame->GetMaximum()*0.85, Form("Chi^2 = %.3f", frame->chiSquare()));
   	TxtChi2->SetTextSize(0.04);
   	TxtChi2->SetTextColor(kRed);
   	frame->addObject(TxtChi2);
	std::cout << " ---> Chi^2 = " << frame->chiSquare() << std::endl;

	TH2 *h_Corr = ResBKGmodel->correlationHist();
	
	// ==== Nbkg ==== //
	RooAbsReal* IntSreg = BKGmodel.createIntegral(M_B0, NormSet(M_B0), Range("SGNregion")); //integrate sgn region
	double Is = IntSreg->getVal();
	std::cout << " - Integral in signal region " << Is << std::endl;
	RooAbsReal* IntBreg = BKGmodel.createIntegral(M_B0, NormSet(M_B0), Range("B0_lSB,B0_rSB"));   //integrate bkg region
	double Ib = IntBreg->getVal();
	std::cout << " - Integral in sidebands region " << Ib << std::endl;

	Nbkg = data->sumEntries() * Is/Ib; 
	std::cout << " - Nbkg in signal region " << Nbkg << std::endl;

	// ==== SAVE ON FILE ==== //
	TString outFileName = Form("./outRoot/PunziOptimization/fit_results/FitResults_B0sidebands_%d_", year_);
	if (BDToutCut_ < 0.) outFileName.Append(Form("m"));
	outFileName.Append(Form("%.0f_Mpipi%.0f", fabs(BDToutCut_)*100., MpipiCut_*1000) + channel_ + ".root");

	TFile* outFitFile = new TFile(outFileName, "RECREATE");
	frame->Write("FitPlot");
	ResBKGmodel->Write("FitResults");
	h_Corr->Write();
	//wsFit->writeToFile("prova_workspace.root");

	outFitFile->Close();
	std::cout << " [OUT] written on " << outFileName << std::endl;
	wsFit->Print();
    return Nbkg;

}//Nbkg_extraction

double MVAoptimizer::Nsgn_extraction(){

	//TString selection = Form("is_signal && M_B0 > %f && M_B0 < %f && BDTout > %f && M_PiPi > %f", B0_SR_low, B0_SR_high, BDToutCut_, MpipiCut_ );
	TString selection = SGNselection + Form(" && M_B0 > %f && M_B0 < %f", B0_SR_low, B0_SR_high); 
    //if(channel_ == "X3872") selection.Append("&& M_X3872 >3.75");
    //lse if(channel_ == "Psi2S") selection.Append("&& M_X3872 < 3.75");
	double Nsgn = inTree_->Draw("M_B0", selection, "goff");
	std::cout << " - Nsgn in signal region " <<  Nsgn << std::endl;

	Nsgn *= LumiNorm_;
	return Nsgn;

}//Nsgn_extraction()

double MVAoptimizer::EFFsgn_extraction(){

	double Ns = Nsgn_extraction();
	TString selection = Form("is_signal && M_B0 > %f && M_B0 < %f", B0_SR_low, B0_SR_high);
    //if(channel_ == "X3872") selection.Append("&& M_X3872 >3.75");
    //else if(channel_ == "Psi2S") selection.Append("&& M_X3872 < 3.75");
	double Ntot = inTree_->Draw("M_B0", selection, "goff");
	Ntot *= LumiNorm_;

	return Ns/Ntot;

}//EFFsgn_extraction()

double MVAoptimizer::PunziSign(double* PSerr){

	double PunziS;
	const double b = 5.0;    // #sigmas corresp. to 5sigma significance level
	const double a = 2.0;    // #sigmas corresp. to CL (90%--> 1.2816) (95% --> 1.6448) (CMStwiki --> 2.)

	double B          = Nbkg_extraction();

	double S          = Nsgn_extraction();
	double Seff       = EFFsgn_extraction();
	double Seff_error = 1./sqrt(S); // error square
	
	double Sign_denom        = b*b + 2.*a*sqrt(B) + b*sqrt(b*b + 4.*a*sqrt(B) + 4.*B );
	double Sign_denom_error  = a/sqrt(B) + b * (a / sqrt(B) + 2. ) / sqrt(b*b + 4.*a*sqrt(B) + 4.*B );	
	Sign_denom_error /= Sign_denom;	

	PunziS = Seff/Sign_denom;
	*PSerr = PunziS*sqrt(Seff_error + Sign_denom_error*Sign_denom_error);

	return PunziS;

}//PunziSign()

int MVAoptimizer::makeSGNvsBKGplot(){

	int Nbins = 45;
	double Mlow = MB_low, Mhigh = MB_high;
	TH1F* h_Data_B0 = new TH1F("Data_B0", "", Nbins, Mlow, Mhigh);
	TH1F* h_SGN_B0  = new TH1F("SGN_B0" , "", Nbins, Mlow, Mhigh); 
	double Mlow_X = 3.75, Mhigh_X = 4.0;
	if (channel_ == "Psi2S") { Mlow_X = 3.60, Mhigh_X = 3.75; }
	TH1F* h_Data_X = new TH1F("Data_X", "", 30, Mlow_X, Mhigh_X);
	TH1F* h_SGN_X  = new TH1F("SGN_X" , "", 30, Mlow_X, Mhigh_X); 
	TH2F* h_Data_B0vsX = new TH2F("Data_B0vsX", "", Nbins, Mlow, Mhigh, 15, Mlow_X, Mhigh_X );
	TH2F* h_SGN_B0vsX  = new TH2F("SGN_B0vsX", "", Nbins, Mlow, Mhigh, 15, Mlow_X, Mhigh_X );

	TString BKG_selection = DATAselection;//Form("!is_signal&& BDTout > %f && M_PiPi > %f",  BDToutCut_, MpipiCut_ );
	if (is_blind_) BKG_selection.Append(Form("&& (M_B0 < %f || M_B0 > %f)",  B0_lSB_high, B0_rSB_low));
	TString SGN_selection = SGNselection;//Form("is_signal && BDTout > %f && M_PiPi > %f",  BDToutCut_, MpipiCut_ );
	inTree_->Draw("M_B0>>Data_B0", BKG_selection);
	inTree_->Draw("M_X3872>>Data_X", BKG_selection+ Form("&& M_B0 > %f && M_B0 < %f", B0_SR_low, B0_SR_high));
	inTree_->Draw("M_X3872:M_B0>>Data_B0vsX", BKG_selection);
	inTree_->Draw("M_B0>>SGN_B0", SGN_selection);
    std::cout << " #sgn @ plot func " << h_SGN_B0->Integral() << std::endl;
	h_SGN_B0->Scale(LumiNorm_);
	inTree_->Draw("M_X3872>>SGN_X", SGN_selection);
	h_SGN_X->Scale(LumiNorm_);
	inTree_->Draw("M_X3872:M_B0>>SGN_B0vsX", SGN_selection);
	h_SGN_B0vsX->Scale(LumiNorm_);
    h_SGN_B0vsX->Add(h_Data_B0vsX);

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


    h_SGN_B0vsX->GetXaxis()->SetTitle("M(B^{0}) [GeV]");
    h_SGN_B0vsX->GetYaxis()->SetTitle("M(J#psi #pi^{+}#pi^{-}) [GeV]");
    h_SGN_B0vsX->GetXaxis()->SetTitleOffset(1.1); h_SGN_B0vsX->GetXaxis()->SetTitleSize(0.04);
	h_SGN_B0vsX->GetYaxis()->SetTitleSize(0.04);
	h_SGN_B0vsX->SetLineWidth(3);
	h_SGN_B0vsX->SetLineColor(kAzure + 1); //h_SGN_B0vsX->SetFillColorAlpha(kAzure + 1, 0.30);
	h_Data_B0vsX->SetLineColor(kRed);
	h_Data_B0vsX->SetLineWidth(3);

	// get FIT
	RooCurve* FitCurve = (RooCurve*)((RooPlot*)wsFit->obj("FitPlot"))->getObject(1);//(RooCurve*)FitPlot->getObject(1);
	FitCurve->SetLineColor(kRed);

	// Legend 
	auto legendB0 = new TLegend(0.53, 0.70,.83,.83);
	legendB0->SetTextSize(0.035);
	legendB0->SetBorderSize(0);

	gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
	gStyle->SetOptStat(0);
	gStyle->SetLineWidth(3);

	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
	h_SGN_B0->GetYaxis()->SetRangeUser(0., 1.4 * std::max(h_SGN_B0->GetMaximum(), h_Data_B0->GetMaximum()));
	h_SGN_B0->Draw("HIST");
	legendB0->AddEntry(h_SGN_B0, "SIGNAL (MC)");
	legendB0->AddEntry(h_Data_B0, Form("DATA %d", year_));
	FitCurve->Draw("SAME");
	h_Data_B0->Draw("PE0 SAME");
	legendB0->AddEntry(FitCurve, "SIDEBANDS-FIT");

	legendB0->Draw();
	//TString outPath = "./outRoot/PunziOptimization/";
	CMSxxx(c1);
	c1->SaveAs(outPath_ + "B0massCUT_" + Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000) + channel_ + ".png");
	c1->SaveAs(outPath_ + "B0massCUT_" + Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000.)+ channel_ + ".pdf");

	h_SGN_X->GetYaxis()->SetRangeUser(0., 1.2 * std::max(h_SGN_X->GetMaximum(), h_Data_X->GetMaximum()));
	h_SGN_X->Draw("HIST");
	h_Data_X->Draw("PE0 SAME");
    //c1 = myRootLib::RatioPlot(h_Data_X,h_SGN_X);
	CMSxxx(c1);
	c1->SaveAs(outPath_ + "JpsiPiPimassCUT_" + Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000) + channel_ + ".png");
	c1->SaveAs(outPath_ + "JpsiPiPimassCUT_" + Form("X%.0f_Mpipi%.0f_", fabs(BDToutCut_)*100., MpipiCut_*1000.)+ channel_ + ".pdf");


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
	legend->AddEntry(h_CorrD, Form("DATA %d Corr = %.3f", year_, h_CorrD->GetCorrelationFactor()));
	h_CorrS->Draw("BOX ");
	h_CorrD->Draw("BOX SAME");
	h_CorrS->Draw("BOX SAME");
    Mpipi_th.Draw("same"); BDTout_th.Draw("same");
	legend->Draw();
    CMSxxx(c1);

	c1->SaveAs(outPath_ + Form("BDToutVSMpipi_UL%d_", year_) + channel_+ ".png"); 
	c1->SaveAs(outPath_ + Form("BDToutVSMpipi_UL%d_", year_) + channel_+ ".pdf"); 

	h_Mpipi_B->Draw("HIST"); h_Mpipi_S->Draw("HIST SAME");
	gPad->SetLeftMargin(0.15);
	legend->Draw();
	c1->SaveAs(outPath_ + Form("Mpipi_dataVSmc_UL%d_", year_) + channel_+ ".png"); 
	c1->SaveAs(outPath_ + Form("Mpipi_dataVSmc_UL%d_", year_) + channel_+ ".pdf"); 

	h_BDTx_B->Draw("HIST"); h_BDTx_S->Draw("HIST SAME");
	gPad->SetLeftMargin(0.15);
	legend->Draw();
	c1->SaveAs(outPath_ + Form("BDTout_dataVSmc_UL%d_", year_) + channel_+ ".png"); 
	c1->SaveAs(outPath_ + Form("BDTout_dataVSmc_UL%d_", year_) + channel_+ ".pdf"); 
	

}

void MVAoptimizer::CMSxxx(TCanvas* c){
	c->cd();
	TLatex RunDetails; RunDetails.SetNDC(); 
	RunDetails.SetTextFont(61);
	RunDetails.SetTextAlign(10);
	RunDetails.SetTextSize(0.035);
	RunDetails.DrawLatex(.12, .92, "CMS");
	RunDetails.SetTextFont(52);
	RunDetails.SetTextSize(0.035);
	RunDetails.DrawLatex(.20, .91, "Work in progress");
	RunDetails.SetTextFont(42);
	RunDetails.SetTextSize(0.030);
	RunDetails.DrawLatex(.70, .91, "41 fb^{-1} (13 TeV)");

}
