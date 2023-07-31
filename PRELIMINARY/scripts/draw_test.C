void draw_test(std::string tag = ""){
  
  // Misc settings
  gROOT->SetBatch(kTRUE);
  gStyle->SetLineScalePS(1.75);
  gStyle->SetOptStat(0);

  // Set colors
  Int_t Number = 5;
  Double_t Stops[5] = {0, .25, .50, .75, 1.};
  Double_t Red[5] = {0, 1., 1., 247./255., 61./255.};
  Double_t Green[5] = {162./255., 0.455, 240./255., 121./255., 1.};
  Double_t Blue[5] = {1., 0., 33./255., 125./255., 109./255.};  
  TColor::CreateGradientColorTable(Number, Stops, Red, Green, Blue, 255, 0.5);

  // gStyle->SetPalette(kLightTemperature, 0, 0.5);

  // Opening files
  std::string outpath = "~/analysisB0ToX3872K0s/PRELIMINARY/my_plots/noemi/";

  TFile* f;
  Bool_t is_X = tag == "X";
  if(!is_X) f = TFile::Open("root://xrootd-cms.infn.it///store/user/crovelli/Run32022__BuToPsi2SKs_Psi2SToJPsiPiPi_JPsiToMuMu_RhoToPiPi_TuneCP5_13p6TeV_pythia8-evtgen/BuToPsi2SKs_Psi2SToJPsiPiPi_JPsiToMuMu_RhoToPiPi.root");
  else if(is_X) f = TFile::Open("root://xrootd-cms.infn.it///store/user/crovelli/Run32022__BuToX3872Ks_X3872ToJPsiRho_JPsiToMuMu_RhoToPiPi_TuneCP5_13TeV_pythia8-evtgen/BuToX3872Ks_X3872ToJPsiRho_JPsiToMuMu_RhoToPiPi.root");
  if(is_X) tag += "_";

  std::string title_dataset_append = " for B^{  0} #rightarrow #psi(2S) K^{0}_{s}";
  if(is_X) title_dataset_append = " for B^{  0} #rightarrow X(3872) K^{0}_{s}";

  TTree* evs = (TTree*)f->Get("Events");


  // Number of secondary vertices
  TCanvas* c1 = new TCanvas();
  evs->Draw("nSV >> hnSV", "", "PLC PFC");
  ((TH1F*)gDirectory->Get("hnSV"))->SetTitle(("Number of secondary vertices" + title_dataset_append).c_str());
  ((TH1F*)gDirectory->Get("hnSV"))->GetXaxis()->SetTitle("# SV");
  ((TH1F*)gDirectory->Get("hnSV"))->GetYaxis()->SetTitle("Events");
  c1->SaveAs((outpath + tag + "n_secondaryVertices.pdf").c_str());

  // Particle multiplicities
  TCanvas* c2 = new TCanvas();
  evs->Draw("nB0 >> hnB0", "", "PLC PFC");
  evs->Draw("nK0s >> hnK0s", "", "PLC PFC SAMES");
  evs->Draw("npipi >> hnpipi", "", "PLC PFC SAMES");
  evs->Draw("nJPsiToMuMu >> hnJPsiToMuMu", "", "PLC PFC SAMES");
  evs->Draw("nMuon >> hnMuon", "", "PLC PFC SAMES");

  gPad->SetLeftMargin(0.125);
  gPad->SetRightMargin(0.075);

  TH1F* hnb0 = (TH1F*)gDirectory->Get("hnB0");
  hnb0->SetTitle(("Particle multiplicities per event" + title_dataset_append).c_str());
  hnb0->SetMaximum(10000);
  hnb0->GetXaxis()->SetRangeUser(0, 35);
  hnb0->GetXaxis()->SetTitle("Number of particles");
  hnb0->GetYaxis()->SetTitle("Events");

  auto leg2 = new TLegend(.45, .55, .925, .9);
  leg2->SetTextSize(.05);
  leg2->AddEntry("hnB0", TString::Format("B^{ 0}: %.2f #pm %.2f", ((TH1F*)gDirectory->Get("hnB0"))->GetMean(), ((TH1F*)gDirectory->Get("hnB0"))->GetStdDev()), "f");
  leg2->AddEntry("hnK0s", TString::Format("K^{0}_{s}: %.2f #pm %.2f", ((TH1F*)gDirectory->Get("hnK0s"))->GetMean(), ((TH1F*)gDirectory->Get("hnK0s"))->GetStdDev()), "f");
  leg2->AddEntry("hnpipi", TString::Format("#pi#pi: %.2f #pm %.2f", ((TH1F*)gDirectory->Get("hnpipi"))->GetMean(), ((TH1F*)gDirectory->Get("hnpipi"))->GetStdDev()), "f");
  leg2->AddEntry("hnJPsiToMuMu", TString::Format("J/#psi #rightarrow #mu#mu: %.2f #pm %.2f", ((TH1F*)gDirectory->Get("hnJPsiToMuMu"))->GetMean(), ((TH1F*)gDirectory->Get("hnJPsiToMuMu"))->GetStdDev()), "f");
  leg2->AddEntry("hnMuon", TString::Format("#mu: %.2f #pm %.2f", ((TH1F*)gDirectory->Get("hnMuon"))->GetMean(), ((TH1F*)gDirectory->Get("hnMuon"))->GetStdDev()), "f");
  leg2->Draw();
  // c2->BuildLegend();

  c2->SaveAs((outpath + tag + "particle_multiplicities.pdf").c_str());

  // // Mixed kinematics
  // TCanvas* c3 = new TCanvas();
  // c3->Divide(3,4);

  // c3->cd(1);
  // evs->Draw("B0_finalFit_mu1_pt", "", "PLC");
  // evs->Draw("B0_finalFit_mu2_pt", "", "PLC SAMES");
  // gPad->BuildLegend();

  // c3->cd(2);
  // evs->Draw("B0_finalFit_mu1_eta", "", "PLC");
  // evs->Draw("B0_finalFit_mu2_eta", "", "PLC SAMES");
  // gPad->BuildLegend();

  // c3->cd(3);
  // evs->Draw("B0_finalFit_pi1_pt", "", "PLC");
  // evs->Draw("B0_finalFit_pi2_pt", "", "PLC SAMES");
  // gPad->BuildLegend();

  // c3->cd(4);
  // evs->Draw("B0_finalFit_JPsi_mass");

  // c3->cd(5);
  // evs->Draw("B0_finalFit_Rho_mass");

  // c3->cd(6);
  // evs->Draw("B0_finalFit_X_mass");

  // c3->cd(7);
  // evs->Draw("B0_K0s_nmcFitted_mass");

  // c3->cd(8);
  // evs->Draw("B0_K0s_nmcFitted_mass");

  // c3->cd(9);
  // evs->Draw("B0_finalFit_pt");

  // c3->cd(10);
  // evs->Draw("B0_finalFit_mass");

  // c3->SaveAs((outpath + tag + "misc_kinematics.pdf").c_str());
}