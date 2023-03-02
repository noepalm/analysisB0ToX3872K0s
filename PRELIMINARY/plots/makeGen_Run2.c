{
	gStyle->SetLineWidth(3);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

    
	// HISTO TRANSVERSE MOMENTUM
	compare_years("gen_pt_mu", "p_{T} (#mu)");
	compare_years("gen_pt_pi", "p_{T} (#pi)");
	compare_years("gen_pt_JPsi", "p_{T} (#mu^{+} #mu^{-})");
	compare_years("gen_pt_Rho", "p_{T} (#pi^{+} #pi^{-})");
	compare_years("gen_pt_X", "p_{T} (#mu^{+} #mu^{-}#pi^{+} #pi^{-})");
	compare_years("gen_pt_Ks", "p_{T} (K_{s}^{0}");
	compare_years("gen_pt_B", "p_{T} (B^{0}");


	compare_years("gen_eta_mu", "#eta (#mu)");
	compare_years("gen_eta_pi", "#eta (#pi)");
	compare_years("gen_eta_JPsi", "#eta (#mu^{+} #mu^{-})");
	compare_years("gen_eta_Rho", "#eta (#pi^{+} #pi^{-})");
	compare_years("gen_eta_X", "#eta (#mu^{+} #mu^{-}#pi^{+} #pi^{-})");
	compare_years("gen_eta_Ks", "#eta (K_{s}^{0}");
	compare_years("gen_eta_B", "#eta (B^{0}");


}
