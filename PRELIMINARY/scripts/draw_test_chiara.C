void draw_test_chiara(std::string tag = ""){
  
  // Misc settings
  gROOT->SetBatch(kTRUE);
  gStyle->SetLineScalePS(1.75);
  gStyle->SetOptStat(0);

  // Set colors
  Int_t Number = 2;
  Double_t Stops[5] = {0, .25, .50, .75, 1.};
  Double_t Red[5] = {0, 1., 1., 247./255., 61./255.};
  Double_t Green[5] = {162./255., 0.455, 240./255., 121./255., 1.};
  Double_t Blue[5] = {1., 0., 33./255., 125./255., 109./255.};  
  TColor::CreateGradientColorTable(Number, Stops, Red, Green, Blue, 255, 0.5);

  std::string outpath = "~/analysisB0ToX3872K0s/PRELIMINARY/my_plots/chiara/";

  TFile *f, *f_HLT;
  if(tag == "") {
    f     = TFile::Open("~/analysisB0ToX3872K0s/PRELIMINARY/outRoot/RecoDecay_Psi2S_Run3.root"); 
    f_HLT = TFile::Open("~/analysisB0ToX3872K0s/PRELIMINARY/outRoot/RecoDecay_Psi2S_HLT_emulation_check.root");
  }
  else if(tag == "X") {
    f = TFile::Open("~/analysisB0ToX3872K0s/PRELIMINARY/outRoot/RecoDecay_X3872_Run3.root"); 
    f_HLT = TFile::Open("~/analysisB0ToX3872K0s/PRELIMINARY/outRoot/RecoDecay_X3872_HLT_emulation_check.root"); 
  }

  //------ HISTOGRAMS ------//
  auto mumu_m_mc = (TH1F*)f->Get("MuMu_M_MC");
  auto mumu_m_fk = (TH1F*)f->Get("MuMu_M_Fk");
  auto mumu_m_prefit = (TH1F*)f->Get("MuMu_M_prefit");
  
  auto pipi_m_mc = (TH1F*)f->Get("PiPi_M_MC");
  auto pipi_m_fk = (TH1F*)f->Get("PiPi_M_Fk");
  auto pipi_m_prefit = (TH1F*)f->Get("PiPi_M_prefit");
    
  auto k0s_m_mc = (TH1F*)f->Get("K0s_M_MC");
  auto k0s_m_fk = (TH1F*)f->Get("K0s_M_Fk");
  auto k0s_m_prefit = (TH1F*)f->Get("K0s_M_prefit");
  
  auto b0_m_mc = (TH1F*)f->Get("B0_M_MC");
  auto b0_m_fk = (TH1F*)f->Get("B0_M_Fk");
  auto b0_m_prefit = (TH1F*)f->Get("B0_M_prefit");
  
  auto x_m_mc = (TH1F*)f->Get("X3872_M_MC");
  auto x_m_fk = (TH1F*)f->Get("X3872_M_Fk");
  auto x_m_prefit = (TH1F*)f->Get("X3872_M_prefit");

  auto pi_pt_mc = (TH1F*)f->Get("Pi1_pT_MC");
  auto pi_pt_fk = (TH1F*)f->Get("Pi1_pT_Fk"); 

  auto pipi_pt_mc = (TH1F*)f->Get("PiPi_pT_MC");
  auto pipi_pt_fk = (TH1F*)f->Get("PiPi_pT_Fk"); 

  auto b0_pt_mc = (TH1F*)f->Get("B0_pT_MC");
  auto b0_pt_fk = (TH1F*)f->Get("B0_pT_Fk"); 

  // trigger emulation
  auto h_trigger_fired = (TH1F*)f_HLT->Get("trigger_fired");
  auto h_trigger_fired_emulated = (TH1F*)f_HLT->Get("trigger_fired_emulated");
  


  // --------- PLOTTING --------- //
  Double_t pad_width = 300;
  Int_t nrows = 3, ncols = 5;
  TCanvas* c1 = new TCanvas("c1", "c1", pad_width*ncols, pad_width*nrows);
  c1->Divide(ncols, nrows);
  TLegend* leg;

  c1->cd(1);
  mumu_m_mc->SetTitle("Mass of #mu#mu");
  mumu_m_mc->GetXaxis()->SetTitle("m(#mu#mu) [GeV]");
  mumu_m_mc->Draw("plc pfc");
  mumu_m_fk->Draw("plc pfc SAMES");

  leg= new TLegend();
  leg->SetTextSize(.035);
  leg->AddEntry(mumu_m_mc, "MC-matched", "f");
  leg->AddEntry(mumu_m_fk, "Fake", "f");
  leg->Draw();
  
  c1->cd(2);
  pipi_m_mc->SetTitle("Mass of #pi#pi");
  pipi_m_mc->GetXaxis()->SetTitle("m(#pi#pi) [GeV]");
  pipi_m_mc->Draw("plc pfc");
  pipi_m_fk->Draw("plc pfc SAMES");
  leg= new TLegend();
  leg->SetTextSize(.035);
  leg->AddEntry(pipi_m_mc, "MC-matched", "f");
  leg->AddEntry(pipi_m_fk, "Fake", "f");
  leg->Draw();
  
  c1->cd(3);
  k0s_m_mc->SetTitle("Mass of K^{0}_{s}");
  k0s_m_mc->GetXaxis()->SetTitle("m(K^{0}_{s}) [GeV]");
  k0s_m_mc->Draw("plc pfc");
  k0s_m_fk->Draw("plc pfc SAMES");
  leg= new TLegend();
  leg->SetTextSize(.035);
  leg->AddEntry(k0s_m_mc, "MC-matched", "f");
  leg->AddEntry(k0s_m_fk, "Fake", "f");
  leg->Draw();

  c1->cd(4);
  std::string x_title = tag == "" ? "m(#psi(2S))" : "m(X(3872))";
  x_m_mc->SetTitle(("Mass of " + x_title.substr(2, x_title.size()-3)).c_str());
  x_m_mc->GetXaxis()->SetTitle((x_title + " [GeV]").c_str());
  x_m_mc->Draw("plc pfc");
  x_m_fk->Draw("plc pfc SAMES");
  leg= new TLegend();
  leg->SetTextSize(.035);
  leg->AddEntry(x_m_mc, "MC-matched", "f");
  leg->AddEntry(x_m_fk, "Fake", "f");
  leg->Draw();

  c1->cd(5);
  b0_m_mc->SetTitle("Mass of B^{ 0}");
  b0_m_mc->GetXaxis()->SetTitle("m(B^{0}) [GeV]");
  b0_m_mc->Draw("plc pfc");
  b0_m_fk->Draw("plc pfc SAMES");
  leg= new TLegend();
  leg->SetTextSize(.035);
  leg->AddEntry(b0_m_mc, "MC-matched", "f");
  leg->AddEntry(b0_m_fk, "Fake", "f");
  leg->Draw();

  c1->cd(6);
  TH1F* mumu_m_mc_copy = (TH1F*)mumu_m_mc->Clone();
  mumu_m_mc_copy->SetTitle((mumu_m_mc->GetTitle() + std::string(", MC-matched")).c_str());
  mumu_m_mc_copy->Draw("plc pfc");
  mumu_m_prefit->Draw("plc pfc SAMES");
  leg= new TLegend();
  leg->SetTextSize(.035);
  leg->AddEntry(mumu_m_mc_copy, "Post-fit", "f");
  leg->AddEntry(mumu_m_prefit, "Pre-fit", "f");
  leg->Draw();
  
  c1->cd(7);
  TH1F* pipi_m_mc_copy = (TH1F*)pipi_m_mc->Clone();
  pipi_m_mc_copy->SetTitle((pipi_m_mc->GetTitle() + std::string(", MC-match")).c_str());
  pipi_m_mc_copy->Draw("plc pfc");
  pipi_m_prefit->Draw("plc pfc SAMES");
  leg= new TLegend();
  leg->SetTextSize(.035);
  leg->AddEntry(pipi_m_mc_copy, "Post-fit", "f");
  leg->AddEntry(pipi_m_prefit, "Pre-fit", "f");
  leg->Draw();
  
  c1->cd(8);
  TH1F* k0s_m_mc_copy = (TH1F*)k0s_m_mc->Clone();
  k0s_m_mc_copy->SetTitle((k0s_m_mc->GetTitle() + std::string(", MC-matched")).c_str());
  k0s_m_mc_copy->Draw("plc pfc");
  k0s_m_prefit->Draw("plc pfc SAMES");
  leg= new TLegend();
  leg->SetTextSize(.035);
  leg->AddEntry(k0s_m_mc_copy, "Post-fit", "f");
  leg->AddEntry(k0s_m_prefit, "Pre-fit", "f");
  leg->Draw();

  c1->cd(9);
  TH1F* x_m_mc_copy = (TH1F*)x_m_mc->Clone();
  x_m_mc_copy->SetTitle((x_m_mc->GetTitle() + std::string(", MC-matched")).c_str());
  x_m_mc_copy->Draw("plc pfc");
  x_m_prefit->Draw("plc pfc SAMES");
  leg= new TLegend();
  leg->SetTextSize(.035);
  leg->AddEntry(x_m_mc_copy, "Post-fit", "f");
  leg->AddEntry(x_m_prefit, "Pre-fit", "f");
  leg->Draw();

  c1->cd(10);
  TH1F* b0_m_mc_copy = (TH1F*)b0_m_mc->Clone();
  b0_m_mc_copy->SetTitle((b0_m_mc->GetTitle() + std::string(", MC-matched")).c_str());
  b0_m_mc_copy->Draw("plc pfc");
  b0_m_prefit->Draw("plc pfc SAMES");
  leg= new TLegend();
  leg->SetTextSize(.035);
  leg->AddEntry(b0_m_mc_copy, "Post-fit", "f");
  leg->AddEntry(b0_m_prefit, "Pre-fit", "f");
  leg->Draw();

  c1->cd(11);
  pipi_pt_mc->SetMaximum(std::max(pipi_pt_mc->GetMaximum(), pipi_pt_fk->GetMaximum())*1.1);
  pipi_pt_mc->SetTitle("Transverse momentum of  #pi#pi");
  pipi_pt_mc->GetXaxis()->SetTitle("p_{T}(#pi#pi)");  
  pipi_pt_mc->Draw("plc pfc");
  pipi_pt_fk->Draw("plc pfc SAMES");
  leg = new TLegend(0.5, 0.7, 0.9, 0.9);
  leg->AddEntry(pipi_pt_mc, "MC-matched", "f");
  leg->AddEntry(pipi_pt_fk, "Fake", "f");
  leg->Draw();

  c1->cd(12);
  b0_pt_mc->SetMaximum(std::max(b0_pt_mc->GetMaximum(), b0_pt_fk->GetMaximum())*1.1);
  b0_pt_mc->SetTitle("Transverse momentum of B^{  0}");
  b0_pt_mc->GetXaxis()->SetTitle("p_{T}(B^{0})");
  b0_pt_mc->Draw("plc pfc");
  b0_pt_fk->Draw("plc pfc SAMES");
  leg = new TLegend(0.5, 0.7, 0.9, 0.9);
  leg->AddEntry(b0_pt_mc, "MC-matched", "f");
  leg->AddEntry(b0_pt_fk, "Fake", "f");
  leg->Draw();

  if(tag != "") tag += "_";
  c1->SaveAs((outpath + tag + "misc_kinematics.pdf").c_str());
  
  // trigger emulation check
  auto c2 = new TCanvas();
  gPad->SetLogy();

  h_trigger_fired->SetMaximum(std::max(h_trigger_fired->GetMaximum(), h_trigger_fired_emulated->GetMaximum())*1.1);
  h_trigger_fired->SetTitle("HLT_DoubleMu4_3_LowMass trigger bit distribution, emulation on MC-matched reco objects");
  h_trigger_fired->GetXaxis()->SetTitle("Trigger bit");
  h_trigger_fired->Draw("PLC PFC");
  h_trigger_fired_emulated->Draw("PLC PFC SAMES");

  // change title size
  gPad->Update(); //force painting to generate title primitive
  TPaveText* pt = (TPaveText*)(c2->GetPrimitive("title"));
  pt->SetTextSize(0.035);
  c2->Modified();

  leg = new TLegend(0.13, 0.6, 0.48, .85);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_trigger_fired, TString::Format("No emulation: %.3f  #pm %.3f", h_trigger_fired->GetMean(), h_trigger_fired->GetStdDev()), "f");
  leg->AddEntry(h_trigger_fired_emulated, TString::Format("Post emulation: %.3f  #pm %.3f", h_trigger_fired_emulated->GetMean(), h_trigger_fired_emulated->GetStdDev()), "f");
  leg->Draw();

  c2->SaveAs((outpath + tag + "trigger_emulation_check.pdf").c_str());

}