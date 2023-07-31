void draw_comparison_plots_m(){
    // misc settings
    gStyle->SetOptStat(0);

    // setup colors
    Int_t Number = 5;
    Double_t Stops[5] = {0, .25, .50, .75, 1.};
    Double_t Red[5] = {0, 1., 1., 247./255., 61./255.};
    Double_t Green[5] = {162./255., 0.455, 240./255., 121./255., 1.};
    Double_t Blue[5] = {1., 0., 33./255., 125./255., 109./255.};  
    TColor::CreateGradientColorTable(Number, Stops, Red, Green, Blue, 255, 0.3);

    // Output folder
    std::string outpath = "/eos/home-n/npalmeri/www/Analysis/preliminary_plots/";

    // Retrive MC files
    TFile* mc_X3872 = TFile::Open("~/analysisB0ToX3872K0s/PRELIMINARY/outRoot/RecoDecay_X3872_Run3.root");
    TFile* mc_Psi2S = TFile::Open("~/analysisB0ToX3872K0s/PRELIMINARY/outRoot/RecoDecay_Psi2S_Run3.root");

    // Retrieve MC trees
    TTree* mc_tree_X3872 = (TTree*)mc_X3872->Get("HLTemulation");
    TTree* mc_tree_Psi2S = (TTree*)mc_Psi2S->Get("HLTemulation");

    // Setup TCanvas for invariant mass plots
    int nrows = 2, ncols = 2;
    TCanvas* c1 = new TCanvas("c1", "c1", 300*nrows, 300*ncols);
    c1->Divide(nrows, ncols);


    // ------ DRAWING MC ------ //

    // Jpsi invariant mass
    c1->cd(1);
    double xlow_jpsi = 3.0, xhigh_jpsi = 3.2;
    TH1F* h_m_jpsi_mc_X3872 = new TH1F("h_m_jpsi_mc_X3872", "h_m_jpsi_mc_X3872", 100, xlow_jpsi, xhigh_jpsi);
    mc_tree_X3872->Draw("M_MuMu >> h_m_jpsi_mc_X3872", "", "PFC PLC HIST");
    h_m_jpsi_mc_X3872->Scale(1./h_m_jpsi_mc_X3872->Integral());
    h_m_jpsi_mc_X3872->SetTitle("m(J/#psi)");
    h_m_jpsi_mc_X3872->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
    h_m_jpsi_mc_X3872->GetYaxis()->SetTitle("dN/dm [GeV^{ -1}]");
    gPad->SetLeftMargin(0.15);

    auto leg1 = new TLegend(0.65, 0.45, 0.85, 0.85);
    leg1->AddEntry(h_m_jpsi_mc_X3872, "MC", "f");

    // Rho = PiPi invariant mass
    c1->cd(2);
    double xlow_rho = 0., xhigh_rho = 1.5;
    TH1F* h_m_rho_mc_X3872 = new TH1F("h_m_rho_mc_X3872", "h_m_rho_mc_X3872", 100, xlow_rho, xhigh_rho);
    mc_tree_X3872->Draw("M_PiPi >> h_m_rho_mc_X3872", "", "PFC PLC HIST");
    h_m_rho_mc_X3872->Scale(1./h_m_rho_mc_X3872->Integral());
    h_m_rho_mc_X3872->SetTitle("m(#pi#pi)");
    h_m_rho_mc_X3872->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
    h_m_rho_mc_X3872->GetYaxis()->SetTitle("dN/dm [GeV^{ -1}]");
    gPad->SetLeftMargin(0.15);

    TH1F* h_m_rho_mc_Psi2S = new TH1F("h_m_rho_mc_Psi2S", "h_m_rho_mc_Psi2S", 100, xlow_rho, xhigh_rho);
    mc_tree_Psi2S->Draw("M_PiPi >> h_m_rho_mc_Psi2S", "", "PFC PLC HIST SAME");
    h_m_rho_mc_Psi2S->Scale(1./h_m_rho_mc_Psi2S->Integral());
    
    auto leg2 = new TLegend(0.6, 0.4, 0.85, 0.85);
    leg2->AddEntry(h_m_rho_mc_X3872, "MC X(3872)", "f");
    leg2->AddEntry(h_m_rho_mc_Psi2S, "MC #psi(2S)", "f");

    // K0s invariant mass
    c1->cd(3);
    double xlow_k0s = 0.4, xhigh_k0s = 0.6;
    TH1F* h_m_k0s_mc_X3872 = new TH1F("h_m_k0s_mc_X3872", "h_m_k0s_mc_X3872", 100, xlow_k0s, xhigh_k0s);
    mc_tree_X3872->Draw("M_K0s >> h_m_k0s_mc_X3872", "", "PFC PLC HIST");
    h_m_k0s_mc_X3872->Scale(1./h_m_k0s_mc_X3872->Integral());
    h_m_k0s_mc_X3872->SetTitle("m(K_{S}^{0})");
    h_m_k0s_mc_X3872->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
    h_m_k0s_mc_X3872->GetYaxis()->SetTitle("dN/dm [GeV^{ -1}]");
    gPad->SetLeftMargin(0.15);

    auto leg3 = new TLegend(0.6, 0.4, 0.85, 0.85);
    leg3->AddEntry(h_m_k0s_mc_X3872, "MC", "f");

    // Psi2S invariant mass
    c1->cd(4);
    double xlow_x3872 = 3.4, xhigh_x3872 = 3.75;
    TH1F* h_m_x3872_mc_Psi2S = new TH1F("h_m_x3872_mc_Psi2S", "h_m_x3872_mc_Psi2S", 100, xlow_x3872, xhigh_x3872);
    mc_tree_Psi2S->Draw("M_X3872 >> h_m_x3872_mc_Psi2S", "", "PFC PLC HIST");
    h_m_x3872_mc_Psi2S->Scale(1./h_m_x3872_mc_Psi2S->Integral());
    h_m_x3872_mc_Psi2S->SetTitle("m(J/#psi #pi #pi)");
    h_m_x3872_mc_Psi2S->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
    h_m_x3872_mc_Psi2S->GetYaxis()->SetTitle("dN/dm [GeV^{ -1}]");
    gPad->SetLeftMargin(0.15);

    auto leg4 = new TLegend(0.2, 0.4, 0.5, 0.85);
    leg4->AddEntry(h_m_x3872_mc_Psi2S, "MC", "f");

    // ------ DRAWING DATA ------ //
    std::vector<std::string> eras = {"D", "E", "F", "G"};
    for(auto era : eras){
        // retrieve data file for corresponding era
        TFile* data_file = TFile::Open(Form("/eos/home-n/npalmeri/B0toX3872K0s/data/Run3_2022%s_blind.root", era.c_str()));

        // retrieve data tree
        TTree* data_tree = (TTree*)data_file->Get("HLTemulation");

        // Jpsi invariant mass
        c1->cd(1);
        TH1F* h_m_jpsi_data = new TH1F("h_m_jpsi_data", "h_m_jpsi_data", 100, xlow_jpsi, xhigh_jpsi);
        data_tree->Draw("M_mumu >> h_m_jpsi_data", "", "PFC PLC HIST SAMES"); //notice different case in variable name
        h_m_jpsi_data->Scale(1./h_m_jpsi_data->Integral());
        leg1->AddEntry(h_m_jpsi_data, Form("Era %s", era.c_str()), "f");

        // Rho = PiPi invariant mass
        c1->cd(2);
        TH1F* h_m_rho_data = new TH1F("h_m_rho_data", "h_m_rho_data", 100, xlow_rho, xhigh_rho);
        data_tree->Draw("M_PiPi >> h_m_rho_data", "", "PFC PLC HIST SAMES");
        h_m_rho_data->Scale(1./h_m_rho_data->Integral());
        leg2->AddEntry(h_m_rho_data, Form("Era %s", era.c_str()), "f");

        // K0s invariant mass
        c1->cd(3);
        TH1F* h_m_k0s_data = new TH1F("h_m_k0s_data", "h_m_k0s_data", 100, xlow_k0s, xhigh_k0s);
        data_tree->Draw("M_K0s >> h_m_k0s_data", "", "PFC PLC HIST SAMES");
        h_m_k0s_data->Scale(1./h_m_k0s_data->Integral());
        leg3->AddEntry(h_m_k0s_data, Form("Era %s", era.c_str()), "f");

        // Psi2S invariant mass
        c1->cd(4);
        TH1F* h_m_x3872_data = new TH1F("h_m_x3872_data", "h_m_x3872_data", 100, xlow_x3872, xhigh_x3872);
        data_tree->Draw("M_X3872 >> h_m_x3872_data", "", "PFC PLC HIST SAMES");
        h_m_x3872_data->Scale(1./h_m_x3872_data->Integral());
        leg4->AddEntry(h_m_x3872_data, Form("Era %s", era.c_str()), "f");
    }

    // Build legend in all subcanvas
    std::vector<TLegend*> legs = {leg1, leg2, leg3, leg4};
    for(int i = 1; i < 5; i++){
        c1->cd(i);
        legs[i-1]->Draw();
    }

    // Save canvas
    c1->SaveAs(Form("%s/compare_mc_data.png", outpath.c_str()));
    c1->SaveAs(Form("%s/compare_mc_data.pdf", outpath.c_str()));
}