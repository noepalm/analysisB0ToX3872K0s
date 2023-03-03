#include <iostream>

#include <TFile.h>
#include <TString.h>

#include <TH1.h>
#include <TH1F.h>
#include <TTree.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>

using namespace std;
//https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html
//This macro contains following functions:
// void draw_single_histogram(const TString histo_name)
// void draw_two_histograms()
// TGraph* makeROCcurve(TH1F* sigHist, TH1F* bkgHist){

TString inRootFile_  = "../outRoot/RecoDecay_X3872_UL17_X3872.root"; 
TString outPath_     = "/eos/user/c/cbasile/www/B0toX3872K0s/HLT-emulation/";


void SetInputFile(const TString& inFile = ""){

  inRootFile_ = inFile;

}//SetInputFile()

void SetOutputFile(const TString& outPath = ""){

  outPath_ = outPath;

}//SetOutputFile()

TFile* open_file(){
    TFile* input_file = new TFile(inRootFile_);

    if ( !input_file->IsOpen() ) {
       std::cout << "ERROR IN OPENING FILE "<< inRootFile_ << std::endl;
       exit(-1);
    }

    return input_file;
}


Color_t PtlColorMap(const TString& particle){

  std::map <TString , Color_t> PtlColor{};
  PtlColor["mu"] = kPink + 5;
  
  PtlColor["muL"] = kAzure + 5;
  PtlColor["muSL"] = kPink + 5;
  PtlColor["pi"] = kBlue - 4;
  PtlColor["PiPi"] = kGreen + 1;
  PtlColor["PiPi_Fk"] = kYellow - 7; 
  PtlColor["JPsi"] = kRed;
  PtlColor["Rho"] = kOrange + 8;
  PtlColor["K0s"] = kGreen;
  PtlColor["K0s_Fk"] = kPink - 6;
  PtlColor["X"] = kBlue - 9;
  PtlColor["X_Fk"] = kMagenta - 6;
  PtlColor["Psi"] = kBlue - 9;
  PtlColor["Psi_Fk"] = kMagenta - 6;
  PtlColor["B0"] = kAzure + 1;
  PtlColor["B0"] = kRed;

  PtlColor["MCmatch"] = kAzure +1;
  PtlColor["MCfake"]  = kRed;
  PtlColor["mu_MC"]   = kBlue;
  PtlColor["mu_Fk"]   = kPink + 5;
  PtlColor["pi_MC"]   = kGreen + 1;
  PtlColor["pi_Fk"]   = kYellow - 7;
  PtlColor["PV"]   = kMagenta -7;
  PtlColor["BS"]   = kGreen -7;

  PtlColor["16preVFP"] = kRed; 
  PtlColor["16"] = kOrange; 
  PtlColor["17"] = kBlue; 
  PtlColor["18"] = kBlack; 

  return PtlColor[particle];
}


TString CategoryLegend(const TString& category){

  std::map <TString , TString> Leg_entry{};
  Leg_entry["mu"] = "#mu";
  Leg_entry["muL"] = "Leading #mu";
  Leg_entry["muSL"] = "Sub-leading #mu";
  Leg_entry["pi"] = "pions";
  Leg_entry["PiPi"] = "#pi^{+} #pi^{-} MC match";
  Leg_entry["PiPi_Fk"] = "#pi^{+} #pi^{-} fake";
  Leg_entry["JPsi"] = "J/#psi MC match";
  Leg_entry["JPsi_Fk"] = "J/#psi fake";
  Leg_entry["Rho"] = "#rho(770)";
  Leg_entry["K0s"] = "K_{s}^{0} MC match";
  Leg_entry["K0s"] = "K_{s}^{0} fake";
  Leg_entry["X"] = "X(3872) MC match";
  Leg_entry["X_Fk"] = "X(3872) fake";
  Leg_entry["Psi"] = "#psi (2S) MC match";
  Leg_entry["Psi_Fk"] = "#psi (2S) fake";
  Leg_entry["B0"] = "B_{0} ";

  Leg_entry["MCmatch"] = "B^{0} MC matching";
  Leg_entry["MCfake"] = "B^{0} FAKE";
  Leg_entry["mu_MC"] = "#mu MC matching";
  Leg_entry["mu_Fk"] = "#mu FAKE";
  Leg_entry["pi_MC"] = "#pi MC matching";
  Leg_entry["pi_Fk"] = "#pi FAKE";
  Leg_entry["PV"] = "wrt chosen PV";
  Leg_entry["BS"] = "wrt BS";

  Leg_entry["16preVFP"] = "2016 pre VFP"; 
  Leg_entry["16"] = "2016 post VFP"; 
  Leg_entry["17"] = "2017"; 
  Leg_entry["18"] = "2018"; 

  return Leg_entry[category];
}

void histoSetUp(TH1* histo, const TString& category, const TString& x_name, bool fill = true , bool norm = true){

  //AXIS LABEL 
  histo->GetXaxis()->SetTitle(x_name);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetLabelSize(0.035);
  histo->GetYaxis()->SetTitle(Form("Events/ %.3f", histo->GetXaxis()->GetBinWidth(1)));
  histo->GetYaxis()->SetTitleOffset(1.6);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetLabelSize(0.035);

  //WIDTH & COLOR
  histo->SetLineWidth(3);
  histo->SetLineColor(PtlColorMap(category));
  if (fill) histo->SetFillColorAlpha(PtlColorMap(category), 0.3);

  //NORMALIZATION
  if(norm )histo->Scale(1./histo->Integral());
}

TString pngName(TString histo_name){

	TString pngName = outPath_ + histo_name + ".png";
	return pngName;
}
TString pdfName(TString histo_name){

	TString pdfName = outPath_ + histo_name + ".pdf";
	return pdfName;
}

void CMSxxx(TCanvas* c){
	c->cd();
	TLatex RunDetails; RunDetails.SetNDC(); 
	RunDetails.SetTextFont(61);
	RunDetails.SetTextAlign(10);
	RunDetails.SetTextSize(0.03);
	RunDetails.DrawLatex(.10, .91, "CMS");
	RunDetails.SetTextFont(52);
	RunDetails.DrawLatex(.17, .91, "Simulation");
	RunDetails.SetTextFont(42);
	RunDetails.SetTextSize(0.025);
	RunDetails.DrawLatex(.70, .91, "41 fb^{-1} (13 TeV)");

}

int draw_one_histo(const TString& histo_name, const TString& category, const TString& x_name, TString out_name){
    
    TFile* input_file = open_file();

    TH1F* h = (TH1F*)input_file->Get(histo_name);
    if ( !h ){
      std::cout<< "null pointer for histogram named " << histo_name << std::endl;
      exit(-1);
    }
    if (out_name == "") out_name = histo_name;

	  histoSetUp(h, category, x_name);

    auto legend = new TLegend(0.60,0.75,.80,.80);
	  legend->SetBorderSize(0);
	  legend->SetTextSize(0.035);
    legend->AddEntry(h,CategoryLegend(category),"f");

    //STATISTICS
    gStyle->SetOptStat(0);

	 //TEXT

    TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);
	  c1->DrawFrame(0,0,1,1);
    h->Draw("HIST");
	  legend->Draw();
	  gPad->SetLeftMargin(0.13);
	  gPad->SetBottomMargin(0.13);
    c1->Update(); 
    c1->SaveAs(pngName(out_name));
    c1->SaveAs(pdfName(out_name));

    input_file->Close();
    return 0;

}//draw_pT()


int draw_two_histograms(const TString histo1,const TString& category1, const TString histo2, const TString& category2, const TString& x_name, TString out_name, bool fill = true, bool LogY = false){

    TFile* input_file = open_file();
    TH1F* h1 = (TH1F*)input_file->Get(histo1);
    TH1F* h2 = (TH1F*)input_file->Get(histo2);
    if ( !h1 ){
      std::cout<< "null pointer for histogram named " << histo1 << std::endl;
      exit(-1);
    }    
    if ( !h2 ){
      std::cout<< "null pointer for histogram named " << histo2 << std::endl;
      exit(-1);
    } 
    if (out_name == "") out_name = histo1;

    histoSetUp(h1, category1, x_name);
    histoSetUp(h2, category2, x_name);
	  h1->GetYaxis()->SetTitle(Form("1/N dN/%.2f GeV", h1->GetXaxis()->GetBinWidth(1)));

    //STATISTICS
    gStyle->SetOptStat(0);
    
    //SETMAXIMUM                                                                                                                                                                  
    double M1 = h1->GetBinContent(h1->GetMaximumBin());
    double M2 = h2->GetBinContent(h2->GetMaximumBin());
    if (M1 > M2){ h1->SetMaximum(1.2*M1);
    }else {h1->SetMaximum(1.2*M2);}
    
    //LEGEND
    auto legend = new TLegend(0.50,0.70,.80,.80);
	  legend->SetBorderSize(0);
	  legend->SetTextSize(0.035);
    legend->AddEntry(h1, CategoryLegend(category1) ,"f");
    legend->AddEntry(h2, CategoryLegend(category2) ,"f");
    
    TString png_name = pngName(out_name);
    TString pdf_name = pdfName(out_name);
    TCanvas* c1 = new TCanvas("c1","canvas", 1024,1024);

    h1->Draw("HIST");
    h2->Draw("HIST SAME");
	  gPad->SetLeftMargin(0.13);
	  gPad->SetBottomMargin(0.13);
    gPad->RedrawAxis();
    legend->Draw();
	  if (LogY) c1->SetLogy();
    else c1->SetLogy(0);
    c1->SaveAs(png_name);
    c1->SaveAs(pdf_name);

    input_file->Close();

    return 0;
}


int draw_binary_histo(const TString h_MCmatch, const TString MC_category, const TString& title,TString out_name){

	  TFile* input_file = open_file();

	  TH1F* h1 = (TH1F*)input_file->Get(h_MCmatch);
	  //TH1F* h2 = (TH1F*)input_file->Get(h_Fake);

	  if ( !h1 ){
		 std::cout<< "null pointer for histogram named " << h_MCmatch << std::endl;
		 exit(-1);
	  }
	  //if ( !h2 ){
		// std::cout<< "null pointer for histogram named " << h_Fake << std::endl;
		// exit(-1);
	  //}
	  
	  TString category1 = MC_category; //, category2 = Fk_category;
	  //SETUP
	  histoSetUp(h1, category1, title);
	  h1->GetYaxis()->SetRangeUser(0.01, 3.);
	  //histoSetUp(h2, category2, title);
    h1->GetXaxis()->SetBinLabel(1, "FALSE"); h1->GetXaxis()->SetBinLabel(2, "TRUE");

	  //STATISTICS
	  gStyle->SetOptStat(0);
	  gStyle->SetPaintTextFormat("1.4f");
	  

	  //LEGEND
	  auto legend = new TLegend(0.62,0.8,.89,.89);
		legend->SetBorderSize(0);
	  legend->AddEntry(h1, CategoryLegend(category1) ,"f");
	  //legend->AddEntry(h2, CategoryLegend(category2) ,"f");

	  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	  gPad->SetLogy();
	  h1->Draw("HIST TEXT0");
	  //h2->Draw("HIST TEXT0 SAME");

	  gPad->RedrawAxis();
	  legend->Draw();

    if (out_name == "") out_name = h_MCmatch;
	  c1->SaveAs(pngName(out_name));
    c1->SaveAs(pdfName(out_name));
	  input_file->Close();

	  return 0;

	}


void makeROCcurve(std::vector<TString> SGNhistos, std::vector<TString> BKGhistos, const TString out_name){

  TFile* input_file = open_file();
  TH1F* sigHist = new TH1F();
  TH1F* bkgHist = new TH1F();
  
  int Nobservables = SGNhistos.size();
  int nbins;  
  float sig_integral = 0, bkg_integral = 0;
  
  TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
  TGraph* vec_graph[Nobservables]; 
  TMultiGraph *mg = new TMultiGraph();

  //LEGEND
  auto legend = new TLegend(0.50,0.15,.80,.30);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.03);

  for (int j = 0; j < Nobservables; j++){

    sigHist = (TH1F*)input_file->Get(SGNhistos[j]);
    bkgHist = (TH1F*)input_file->Get(BKGhistos[j]);

    nbins = sigHist->GetNbinsX();
    sig_integral = sigHist->Integral(1,nbins);
    bkg_integral = bkgHist->Integral(1,nbins);
    std::cout << "Histo number " << j << std::endl;
    std::cout<<" total int  sig: "<<sig_integral<<" bkg: "<<bkg_integral<<std::endl;
    std::vector<float> sigPoints(nbins);
    std::vector<float> bkgPoints(nbins);
    for ( int i = nbins; i > 0; i-- ) {
      float sig_slice_integral = sigHist->Integral(i,nbins);
      float bkg_slice_integral = bkgHist->Integral(i,nbins);
      sigPoints.push_back(sig_slice_integral/sig_integral);
      bkgPoints.push_back(bkg_slice_integral/bkg_integral);

      std::cout<<i<<" "<<sig_slice_integral<<" "<<sig_slice_integral/sig_integral<<" "<<bkg_slice_integral<<" "<<bkg_slice_integral/bkg_integral<<std::endl;
    }
    
    vec_graph[j] = new TGraph(sigPoints.size(),&bkgPoints[0], &sigPoints[0]);
    vec_graph[j]->SetLineWidth(4);
    vec_graph[j]->SetLineColor(2+j);
    legend->AddEntry(vec_graph[j], SGNhistos[j], "l");
    mg->Add(vec_graph[j]);
    
    std::cout <<"\n -------------------------\n" << std::endl;
  } // on observables
  //g->GetXaxis()->SetTitle("signal efficiency"); g->GetYaxis()->SetTitle("background efficiency");

    c1->cd();
    mg->Draw("AL");
    mg->SetTitle("; background efficiency; signal efficiency");
    legend->Draw();
    c1->SaveAs(pngName(out_name));
    c1->SaveAs(pdfName(out_name));
  
  input_file->Close();

}




void compare_years(const TString& tree_name, const TString branch_name, const TString& selection, const int Nbins = 100, const float xlow = 0., const float xhigh = 100, TString x_name = " ", TString out_name = " "){

    //TFile* file_16preVFP = new TFile("outRoot/.root"); 
    TFile* file_16       = new TFile("outRoot/DataBlind_HLTemulation_prova_2016G.root"); 
    TFile* file_17       = new TFile("outRoot/DataBlind_HLTemulation_prova_2017B.root"); 
    TFile* file_18       = new TFile("outRoot/DataBlind_HLTemulation_prova_2018A.root"); 

    //TTree* h_16preVFP = (TTree*)file_16preVFP->Get(tree_name);
    TTree* t_16= (TTree*)file_16->Get(tree_name);
    TH1F*  h_16= new TH1F("h_16", "", Nbins, xlow, xhigh);
    t_16->Draw(branch_name+">>h_16", selection);
    TTree* t_17= (TTree*)file_17->Get(tree_name);
    TH1F*  h_17= new TH1F("h_17", "", Nbins, xlow, xhigh);
    t_17->Draw(branch_name+">>h_17", selection);
    TTree* t_18= (TTree*)file_18->Get(tree_name);
    TH1F*  h_18= new TH1F("h_18", "", Nbins, xlow, xhigh);
    t_18->Draw(branch_name+">>h_18", selection);

    
    //histoSetUp(TH1* histo, const TString& category, const TString& x_name, bool fill = true , bool norm = true)
    //histoSetUp(h_16preVFP, "16preVFP", x_name, false, true);
    if(x_name == " ") x_name = branch_name;
    histoSetUp(h_16      , "16"      , x_name, false, true);
    histoSetUp(h_17      , "17"      , x_name, false, true);
    histoSetUp(h_18      , "18"      , x_name, false, true);

    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);

    //LEGEND                                                                                                                                                                    
    auto legend = new TLegend(0.15,0.65,.40,.80);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.035);
    //legend->AddEntry(h_16preVFP,CategoryLegend("16preVFP"),"f");
    legend->AddEntry(h_16,CategoryLegend("16"),"f");
    legend->AddEntry(h_17,CategoryLegend("17"),"f");
    legend->AddEntry(h_18,CategoryLegend("18"),"f");

    TCanvas* c = new TCanvas("c", "", 800, 600);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.13);
    //h_16preVFP->Draw("hist");
    //h_16->SetMaximum(1.5*std::max( std::max(h_16preVFP->GetBinContent(h_16preVFP->GetMaximumBin()), h_16->GetBinContent(h_16->GetMaximumBin())) , 
     //                                       std::max(h_17->GetBinContent(h_17->GetMaximumBin()), h_18->GetBinContent(h_18->GetMaximumBin()) ) ) ); 
    h_16->SetMaximum(1.5*std::max( std::max(h_16->GetBinContent(h_16->GetMaximumBin()), h_17->GetBinContent(h_17->GetMaximumBin())) , 
                                            std::max(h_17->GetBinContent(h_17->GetMaximumBin()), h_18->GetBinContent(h_18->GetMaximumBin()) ) ) ); 
    h_16->Draw("hist");
    h_17->Draw("same hist");
    h_18->Draw("same hist");
    legend->Draw();
    if(out_name == " ") out_name = branch_name;
    c->SaveAs("/eos/user/c/cbasile/www/B0toX3872K0s/HLT-emulation/FullRun2_data/"+ out_name + ".png");
    c->SaveAs("/eos/user/c/cbasile/www/B0toX3872K0s/HLT-emulation/FullRun2_data/"+ out_name + ".pdf");


    //file_16preVFP->Close();
    file_16->Close();
    file_17->Close();
    file_18->Close();

}//compare_years()
