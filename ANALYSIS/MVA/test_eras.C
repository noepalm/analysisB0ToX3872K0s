void test_eras(){
    // disable stats box
    gStyle->SetOptStat(0);

    // Set colors
    Int_t Number = 4;
    Double_t Stops[5] = {0, .25, .50, .75, 1.};
    Double_t Red[5] = {0, 1., 1., 247./255., 61./255.};
    Double_t Green[5] = {162./255., 0.455, 240./255., 121./255., 1.};
    Double_t Blue[5] = {1., 0., 33./255., 125./255., 109./255.};  
    TColor::CreateGradientColorTable(Number, Stops, Red, Green, Blue, 255, 0.3);

    std::vector<std::string> eras = {"D", "E", "F", "G"};

    TCanvas* c = new TCanvas("c", "c", 1000, 1000);
    c->Divide(2, 2);
    TLegend* leg_cosalpha = new TLegend(0.5, 0.5, 0.8, 0.8);
    TLegend* leg_lxy = new TLegend(0.6, .6, 0.9, 0.9);

    for(int i = 0; i < 4; i++){
        TFile* f = TFile::Open(Form("/eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022%s_BDTapplication_CV3_X3872.root", eras[i].c_str()));
        TTree* t = (TTree*)f->Get(Form("CVtraining_2022%s", eras[i].c_str()));

        std::string draw_options = !i ? "PLC PFC HIST" : "PLC PFC HIST SAME";

        c->cd(1);
        t->Draw(Form("CosAlpha3DBSz_B0>>h_cosalpha%s(100, -1.1, 1.1)", eras[i].c_str()), "", draw_options.c_str());
        leg_cosalpha->AddEntry(Form("h_cosalpha%s", eras[i].c_str()), Form("%s", eras[i].c_str()), "f");

        TH1F* h_cosalpha = (TH1F*)gDirectory->Get(Form("h_cosalpha%s", eras[i].c_str()));
        h_cosalpha->Scale(1./h_cosalpha->Integral());
        if(!i){
            h_cosalpha->SetMaximum(0.95);
            h_cosalpha->SetTitle("cos(#alpha_{B^{0}}), normalized");
        } 

        c->cd(2);
        t->Draw(Form("LxySignBSz_B0 >> h_lxy%s(100, -0.01, 200)", eras[i].c_str()), "", draw_options.c_str());
        leg_lxy->AddEntry(Form("h_lxy%s", eras[i].c_str()), Form("%s", eras[i].c_str()), "f");

        TH1F* h_lxy = (TH1F*)gDirectory->Get(Form("h_lxy%s", eras[i].c_str()));
        h_lxy->Scale(1./h_lxy->Integral());
        if(!i){
            h_lxy->SetTitle("L_{xy}/#sigma_{L_{xy}}, normalized");
        }

        if(i < 2){
            c->cd(3 + i);
            t->Draw("CosAlpha3DBSz_B0 >> h_cosalpha_signal(100, -1.1, 1.1)", "is_signal == 1", "PLC PFC HIST");
            t->Draw("CosAlpha3DBSz_B0 >> h_cosalpha_background(100, -1.1, 1.1)", "is_signal == 0", "PLC PFC HIST SAME");
            gPad->BuildLegend();
        }

    }

    c->cd(1);
    gPad->SetLogy();
    leg_cosalpha->Draw();

    c->cd(2);
    leg_lxy->Draw();

    c->SaveAs("cosalpha_lxysign_era_comparison.pdf");  

}