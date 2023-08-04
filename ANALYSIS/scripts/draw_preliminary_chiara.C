void draw_preliminary_chiara(){
    
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

    // std::string outpath = "~/analysisB0ToX3872K0s/ANALYSIS/my_plots/";
    std::string outpath = "/eos/home-n/npalmeri/www/Analysis/preliminary_plots/data/";

    TFile *f, *f_HLT;
    f     = TFile::Open("~/analysisB0ToX3872K0s/ANALYSIS/outRoot/Parking_Run3_HLTemulation_blind.root"); 
    f_HLT = TFile::Open("~/analysisB0ToX3872K0s/ANALYSIS/outRoot/Parking_HLT_emulation_check_HLTemulation_blind.root");
 
    TFile* f_eras[4];
    TTree* t_eras[4];
    std::vector<std::string> eras = {"D", "E", "F", "G"};
    for(int i = 0; i < 4; i++) {
        f_eras[i] = TFile::Open(Form("/eos/home-n/npalmeri/B0toX3872K0s/data/Run3_2022%s_blind.root", eras[i].c_str()));
        t_eras[i] = (TTree*)f_eras[i]->Get("HLTemulation");
    }

    TFile* f_mc_x = TFile::Open("~/analysisB0ToX3872K0s/PRELIMINARY/outRoot/RecoDecay_X3872_Run3.root");
    TFile* f_mc_psi2s = TFile::Open("~/analysisB0ToX3872K0s/PRELIMINARY/outRoot/RecoDecay_Psi2S_Run3.root");

    //------ HISTOGRAMS ------//
    auto mumu_m_postfit = (TH1F*)f->Get("MuMu_M_postfit");
    auto mumu_m_prefit = (TH1F*)f->Get("MuMu_M_prefit");
    
    auto pipi_m_postfit = (TH1F*)f->Get("PiPi_M_postfit");
    auto pipi_m_prefit = (TH1F*)f->Get("PiPi_M_prefit");

    auto k0s_m_postfit = (TH1F*)f->Get("K0s_M_postfit");
    auto k0s_m_prefit = (TH1F*)f->Get("K0s_M_prefit");
    
    auto b0_m_postfit = (TH1F*)f->Get("B0_M_postfit");
    auto b0_m_prefit = (TH1F*)f->Get("B0_M_prefit");
    
    auto x_m_postfit = (TH1F*)f->Get("X3872_M_postfit");
    auto x_m_prefit = (TH1F*)f->Get("X3872_M_prefit");

    // trigger emulation
    auto h_trigger_fired = (TH1F*)f_HLT->Get("trigger_fired");
    auto h_trigger_fired_emulated = (TH1F*)f_HLT->Get("trigger_fired_emulated");    


    // --------- PLOTTING --------- //
    Double_t pad_width = 300;
    Int_t nrows = 1, ncols = 5;
    TCanvas* c1 = new TCanvas("c1", "c1", pad_width*ncols, pad_width*nrows);
    c1->Divide(ncols, nrows);
    TLegend* leg;

    c1->cd(1);
    mumu_m_postfit->SetTitle("Mass of #mu#mu");
    mumu_m_postfit->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
    mumu_m_postfit->Draw("plc pfc");
    mumu_m_prefit->Draw("plc pfc SAMES");
    leg= new TLegend();
    leg->SetTextSize(.035);
    leg->AddEntry(mumu_m_postfit, "Post-fit", "f");
    leg->AddEntry(mumu_m_prefit, "Pre-fit", "f");
    leg->Draw();
    
    c1->cd(2);
    pipi_m_postfit->SetTitle("Mass of #pi#pi");
    pipi_m_postfit->GetXaxis()->SetTitle("m_{#pi#pi} [GeV]");
    pipi_m_postfit->Draw("plc pfc");
    pipi_m_prefit->Draw("plc pfc SAMES");
    leg= new TLegend();
    leg->SetTextSize(.035);
    leg->AddEntry(pipi_m_postfit, "Post-fit", "f");
    leg->AddEntry(pipi_m_prefit, "Pre-fit", "f");
    leg->Draw();
    
    c1->cd(3);
    k0s_m_postfit->SetTitle("Mass of K^{0}_{s}");
    k0s_m_postfit->GetXaxis()->SetTitle("m_{K^{0}_{s}} [GeV]");
    k0s_m_postfit->Draw("plc pfc");
    k0s_m_prefit->Draw("plc pfc SAMES");
    leg = new TLegend();
    leg->SetTextSize(.035);
    leg->AddEntry(k0s_m_postfit, "Post-fit", "f");
    leg->AddEntry(k0s_m_prefit, "Pre-fit", "f");
    leg->Draw();

    c1->cd(4);
    x_m_postfit->SetTitle("Mass of X(3872), blinded");
    x_m_postfit->GetXaxis()->SetTitle("m_{X(3872)} [GeV]");
    x_m_postfit->Draw("plc pfc");
    x_m_prefit->Draw("plc pfc SAMES");
    leg= new TLegend();
    leg->SetTextSize(.035);
    leg->AddEntry(x_m_postfit, "Post-fit", "f");
    leg->AddEntry(x_m_prefit, "Pre-fit", "f");
    leg->Draw();

    c1->cd(5);
    b0_m_postfit->SetTitle("Mass of B^{0}, blinded");
    b0_m_postfit->GetXaxis()->SetTitle("m_{B^{0}} [GeV]");
    b0_m_postfit->Draw("plc pfc");
    b0_m_prefit->Draw("plc pfc SAMES");
    leg= new TLegend();
    leg->SetTextSize(.035);
    leg->AddEntry(b0_m_postfit, "Post-fit", "f");
    leg->AddEntry(b0_m_prefit, "Pre-fit", "f");
    leg->Draw();

    c1->SaveAs((outpath + "misc_kinematics.pdf").c_str());
    
    // trigger emulation check
    auto c2 = new TCanvas();
    gPad->SetLogy();

    h_trigger_fired->SetMaximum(std::max(h_trigger_fired->GetMaximum(), h_trigger_fired_emulated->GetMaximum())*1.1);
    h_trigger_fired->SetTitle("HLT_DoubleMu4_3_LowMass trigger bit distribution");
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

    c2->SaveAs((outpath + "trigger_emulation_check.pdf").c_str());

    // Number of secondary vertices
    TH1F* hnSV = (TH1F*)f->Get("hnSV");

    TCanvas* c3 = new TCanvas();

    hnSV->Draw("PLC PFC");
    hnSV->SetTitle("Number of secondary vertices");
    hnSV->GetXaxis()->SetTitle("# SV");
    hnSV->GetYaxis()->SetTitle("Events");

    c3->SaveAs((outpath + "n_secondaryVertices.pdf").c_str());

    // Particle multiplicities
    TH1I* hnB0 = (TH1I*)f->Get("hnB0");
    TH1I* hnK0s = (TH1I*)f->Get("hnK0s");
    TH1I* hnpipi = (TH1I*)f->Get("hnpipi");
    TH1I* hnJPsiToMuMu = (TH1I*)f->Get("hnJPsiToMuMu");
    TH1I* hnMuon = (TH1I*)f->Get("hnMuon");

    TCanvas* c4 = new TCanvas("c4", "", 1000, 500);
    c4->Divide(2,1);
    TColor::CreateGradientColorTable(5, Stops, Red, Green, Blue, 255, 0.5);

    c4->cd(1);
    hnB0->SetTitle("Particle multiplicities per event, low range");
    hnB0->SetMaximum(7500);
    hnB0->GetXaxis()->SetRangeUser(1, 15);
    hnB0->GetXaxis()->SetTitle("Number of particles");
    hnB0->GetYaxis()->SetTitle("Events");
    hnB0->Draw("PLC PFC");

    hnK0s->Draw("PLC PFC SAMES");
    hnpipi->Draw("PLC PFC SAMES");
    hnJPsiToMuMu->Draw("PLC PFC SAMES");
    hnMuon->Draw("PLC PFC SAMES");

    gPad->SetLeftMargin(0.125);
    gPad->SetRightMargin(0.075);

    leg = new TLegend(.45, .55, .925, .9);
    leg->SetTextSize(.05);
    leg->AddEntry(hnB0, TString::Format("B^{0}: %.2f #pm %.2f", hnB0->GetMean(), hnB0->GetStdDev()), "f");
    leg->AddEntry(hnK0s, TString::Format("K^{0}_{s}: %.2f #pm %.2f", hnK0s->GetMean(), hnK0s->GetStdDev()), "f");
    leg->AddEntry(hnpipi, TString::Format("#pi#pi: %.2f #pm %.2f", hnpipi->GetMean(), hnpipi->GetStdDev()), "f");
    leg->AddEntry(hnJPsiToMuMu, TString::Format("J/#psi: %.2f #pm %.2f", hnJPsiToMuMu->GetMean(), hnJPsiToMuMu->GetStdDev()), "f");
    leg->AddEntry(hnMuon, TString::Format("#mu: %.2f #pm %.2f", hnMuon->GetMean(), hnMuon->GetStdDev()), "f");
    leg->Draw();
    
    c4->cd(2);
    hnK0s->SetTitle("Particle multiplicities per event, high range");
    hnK0s->SetMaximum(2000);
    hnK0s->GetXaxis()->SetRangeUser(0, 100);
    hnK0s->GetXaxis()->SetTitle("Number of particles");
    hnK0s->GetYaxis()->SetTitle("Events");
    hnK0s->Draw("PLC PFC");

    hnB0->Draw("PLC PFC SAMES");
    hnpipi->Draw("PLC PFC SAMES");
    hnJPsiToMuMu->Draw("PLC PFC SAMES");
    hnMuon->Draw("PLC PFC SAMES");

    gPad->SetLeftMargin(0.125);
    gPad->SetRightMargin(0.075);

    leg = new TLegend(.45, .55, .925, .9);
    leg->SetTextSize(.05);
    leg->AddEntry(hnB0, TString::Format("B^{0}: %.2f #pm %.2f", hnB0->GetMean(), hnB0->GetStdDev()), "f");
    leg->AddEntry(hnK0s, TString::Format("K^{0}_{s}: %.2f #pm %.2f", hnK0s->GetMean(), hnK0s->GetStdDev()), "f");
    leg->AddEntry(hnpipi, TString::Format("#pi#pi: %.2f #pm %.2f", hnpipi->GetMean(), hnpipi->GetStdDev()), "f");
    leg->AddEntry(hnJPsiToMuMu, TString::Format("J/#psi: %.2f #pm %.2f", hnJPsiToMuMu->GetMean(), hnJPsiToMuMu->GetStdDev()), "f");
    leg->AddEntry(hnMuon, TString::Format("#mu: %.2f #pm %.2f", hnMuon->GetMean(), hnMuon->GetStdDev()), "f");
    leg->Draw();

    c4->SaveAs((outpath + "particle_multiplicities.pdf").c_str());


    // BDT VARIABLES

    Int_t Number_bdt = 6;
    Double_t Stops_bdt[6] = {0, .2, .4, .6, .8, 1.};
    Double_t Red_bdt[6] = {0, 1., 1., 61./255., 0.7, 0.5};
    Double_t Green_bdt[6] = {162./255., 0.455, 240./255., 1., 0.6, 0.6};
    Double_t Blue_bdt[6] = {1., 0., 33./255., 109./255., 0.6, 0.6};  
    TColor::CreateGradientColorTable(Number_bdt, Stops_bdt, Red_bdt, Green_bdt, Blue_bdt, 255, 0.3);


    TCanvas* c5 = new TCanvas("c5", "", 1000, 400);
    c5->Divide(5,2);

    TLegend* legs[10];
    for(int i = 0; i < 10; i++) legs[i] = new TLegend();

    // TTree* t = (TTree*)f->Get("HLTemulation");
    std::vector<std::string> features = {"pTM_B0", "LxySignBSz_B0", "SVprob_B0", "CosAlpha3DBSz_B0", "LxySignSV_K0s", "SVprob_PiPi", "pT_PiPi", "pT_Pi1", "DR_B0Pi1", "D0_Pi1"};
    std::vector<std::string> features_latex = {"p_{T}(B^{0})/M(B^{0})", "L_{xy}/#sigma_{L_{xy}} B^{0}", "SV prob B^{0}", "cos(#alpha_{3D}) B^{0}", "L_{xy}/#sigma_{L_{xy}} K^{0}_{S}", "SV prob #rho", "p_{T}(#pi#pi)/p_{T}(B^{0})", "p_{T}(#pi_{1})/p_{T}(B^{0})", "#DeltaR(B^{0}, #pi_{1})", "d_{0}(#pi_{1})"};
    std::vector<std::string> features_range = {"(100, 1.5, 12)", "(100, 3, 30)", "(100, 0, 1)", "(100, -1, 1)", "(100, 0, 800)", "(100, 0, 1)", "", "", "(100, 0, 1)", "(100, 0, 10)"};

    for(int j = 0; j < eras.size(); j++){

        auto era = eras[j];
        TTree* t = t_eras[j];

        for(int i = 0; i < features.size(); i++){
            c5->cd(i + 1);

            std::string draw_options = !j ? "PLC PFC HIST" : "PLC PFC HIST SAME";
            t->Draw((features[i] + ">>" + era + "_" + features[i] + "_hist" + features_range[i]).c_str(), "", draw_options.c_str());
            TH1F* h = (TH1F*)gDirectory->Get((era + "_" + features[i] + "_hist").c_str());
            h->SetTitle((features_latex[i] + ", normalized").c_str());
            h->GetXaxis()->SetTitle(features_latex[i].c_str());
            h->GetYaxis()->SetTitle("Events");
            h->Scale(1./h->Integral());

            if(i == 8) h->SetMaximum(0.05);
            if(i == 9) h->SetMaximum(0.07);

            legs[i]->AddEntry((era + "_" + features[i] + "_hist").c_str(), ("Era " + era).c_str(), "f");

            if(i == 1 || i == 3) gPad->SetLogy(); //set logy for LxySignBSz_B0 and CosAlpha3DBSz_B0
        }
    }

    //also draw MC
    TFile* mc_files[2] = {f_mc_x, f_mc_psi2s};
    for(int j = 0; j < 2; j++){
        TTree* t = (TTree*)mc_files[j]->Get("HLTemulation");

        std::string name = !j ? "_mc_X" : "_mc_Psi2S";
        std::string legend_name = !j ? "MC X(3872)" : "MC #psi(2S)";

        for(int i = 0; i < features.size(); i++){
            c5->cd(i + 1);

            t->Draw((features[i] + ">>" + features[i] + name + "_hist" + features_range[i]).c_str(), "", "PLC PFC HIST SAME");
            TH1F* h = (TH1F*)gDirectory->Get((features[i] + name + "_hist").c_str());
            h->SetFillColorAlpha(kBlack, 0.3);
            h->Scale(1./h->Integral());

            legs[i]->AddEntry((features[i] + name + "_hist").c_str(), legend_name.c_str(), "f");
        }
    }

    for(int i = 0; i < 10; i++){
        c5->cd(i + 1);
        legs[i]->Draw();
    }
    
    c5->SaveAs((outpath + "BDT_inputs.pdf").c_str());

}