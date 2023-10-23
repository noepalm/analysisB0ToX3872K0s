void test_different_training(){

	TString outpath = "/eos/home-n/npalmeri/www/Analysis/MVA/";

     // disable stats box
    gStyle->SetOptStat(0);
    //set batch
    gROOT->SetBatch(kTRUE);

    // Set colors
    Int_t Number = 4;
    Double_t Stops[5] = {0, .25, .50, .75, 1.};
    Double_t Red[5] = {0, 1., 1., 247./255., 61./255.};
    Double_t Green[5] = {162./255., 0.455, 240./255., 121./255., 1.};
    Double_t Blue[5] = {1., 0., 33./255., 125./255., 109./255.};  
    TColor::CreateGradientColorTable(Number, Stops, Red, Green, Blue, 255, 0.3);

    TFile* fD = TFile::Open("/eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022D_BDTapplication_CV3_X3872.root");
    TFile* fE = TFile::Open("/eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022E_BDTapplication_CV3_X3872.root");
    TFile* fF = TFile::Open("/eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022F_BDTapplication_CV3_X3872.root");
    TFile* fG = TFile::Open("/eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022G_BDTapplication_CV3_X3872.root");
    
    //retrive the trees
    TTree* tD = (TTree*)fD->Get("CVtraining_2022D");
    TTree* tE = (TTree*)fE->Get("CVtraining_2022E");
    TTree* tF = (TTree*)fF->Get("CVtraining_2022F");
    TTree* tG = (TTree*)fG->Get("CVtraining_2022G");

	
	// bdtout canvas
    TCanvas* c = new TCanvas("c", "", 400*4, 400);
    c->Divide(4,1);
	// bdtout vs cosalpha canvas
	TCanvas* c2 = new TCanvas("c2", "", 400*4, 400);
    c2->Divide(4,1);

    //make tree pointer vector
    std::vector<TTree*> trees = {tD, tE, tF, tG};
    std::vector<std::string> eras = {"D", "E", "F", "G"};

    TLegend* leg;

    //draw bdtout for signal vs background for each era
    for(int i = 0; i < 4; i++){

        c->cd(i+1);
        trees[i]->Draw(Form("BDTout >> bdtout_sgn_%d(100, -12, 6)", i), "is_signal == 1", "PLC PFC HIST");
		TH1F* h_sgn = (TH1F*)gPad->GetPrimitive(Form("bdtout_sgn_%d", i));

        h_sgn->Scale(1./h_sgn->Integral());
        h_sgn->SetTitle(Form("BDT output for era %s", eras[i].c_str()));
        h_sgn->SetMaximum(0.04);
        h_sgn->GetXaxis()->SetTitle("BDT output");

        trees[i]->Draw(Form("BDTout >> bdtout_bkg_%d(100, -12, 6)", i), "is_signal == 0 && !(M_B0 > 5.20 && M_B0 < 5.34 && M_X3872 > 3.85 && M_X3872 < 3.89)", "PLC PFC HIST SAME"); //data is already blinded, but just in case
        TH1F* h_bkg = (TH1F*)gPad->GetPrimitive(Form("bdtout_bkg_%d", i));
        h_bkg->Scale(1./h_bkg->Integral());

        // make legend
        leg = new TLegend();
        leg->AddEntry(h_sgn, "Signal", "f");
        leg->AddEntry(h_bkg, "Background", "f");
        leg->Draw();

		c2->cd(i+1);
		trees[i]->Draw(Form("CosAlpha3DBSz_B0:BDTout >> 2dhist_%d(100, -12, 6, 100, -1.1, 1.1)", i), "is_signal == 0 && !(M_B0 > 5.18903 && M_B0 < 5.37143 && M_X3872 > 3.82562 && M_X3872 < 3.91839)", "colz"); //data is already blinded, but just in case
		TH2F* hist2d = (TH2F*)gPad->GetPrimitive(Form("2dhist_%d", i));
        hist2d->SetTitle(Form("BDT output VS cosAlpha for era %s", eras[i].c_str()));
        hist2d->GetXaxis()->SetTitle("BDTout");
        hist2d->GetYaxis()->SetTitle("cosAlpha");
    }

    
    c->SaveAs(outpath + "BDTout_per_era.png");
    
    // change n. of colors
    TColor::CreateGradientColorTable(3, Stops, Red, Green, Blue, 255, 0.3);
    c2->SaveAs(outpath + "BDTout_vs_cosAlpha_per_era.png");

	

}
