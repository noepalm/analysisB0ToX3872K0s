{
	gStyle->SetLineWidth(3);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

    
	// HISTO TRANSVERSE MOMENTUM
	compare_years("MuMu_M_MC", "M(J #psi) GeV");
	compare_years("PiPi_M_MC", "M(#pi^+ #pi^-) GeV");
	compare_years("X3872_M_MC", "M(X3872) GeV");
	compare_years("K0s_M_MC", "K^{0}_s GeV");
	compare_years("B0_M_MC", "B^{0} GeV");
	compare_years("B0_M_MC", "B^{0} GeV");

    
    compare_years("PiPi_svProb_MC", "P_{sv}(#pi_{1} #pi_{2})" );
    compare_years("PiPi_pT_MC", "p_{T}(#pi #pi)/p_{T}(B^{0})");
    compare_years("Pi1_pT_MC", "p_{T}(#pi_{1})/p_{T}(B^{0})" );
    compare_years("Pi1_D0_MC", "DCA(#pi_{1}, BS)" );
    compare_years("Pi1_DRwrtB_MC", "#Delta R(#pi_{1}, B^{0})" );
    compare_years("K0s_LxySign_wrtBvtx_MC", "L_{xy}(K^{0_s}; B^{0}_{vtx})/#sigma_{xy}");
    compare_years("B0_cosAlpha2DwrtBSwithZ_MC","cos #alpha (#vec{p}_{B}, L_{B})_{2D}");
    
    compare_years("B0_LxySign_wrtBSwithZ_MC","L_{xy}( BS_{z}, B^{0}_{vtx})/#sigma_{xy}");
    
    compare_years("B0_svProb_MC", "P_{sv}(B^{0})" );
    compare_years("B0_pT_MC", "p_{T}(B^{0})/M(B^{0})");

}
