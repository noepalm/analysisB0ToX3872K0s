{
    // .L plots/plot_library.cc++
    // void compare_years(const TString& tree_name, const TString branch_name, const Tstring& selection, const int Nbins = 100, const float xlow = 0., const float xhigh = 100, const TString& x_name "", const TString& outname = ""){
    compare_years("HLTemulation", "M_B0", "M_X3872 < 3.75", 30, 5., 5.6, "M(B^{0}) GeV", "B0_mass_Psi");
    compare_years("HLTemulation", "M_B0", "M_X3872 > 3.75", 30, 5., 5.6, "M(B^{0}) GeV", "B0_mass_X");

    compare_years("HLTemulation", "M_mumu", "M_X3872 < 3.75", 30, 2.95, 3.25, "M(#mu^{+}#mu^{-}) GeV", "MuMu_mass_Psi");
    compare_years("HLTemulation", "M_mumu", "M_X3872 > 3.75", 30, 2.95, 3.25, "M(#mu^{+}#mu^{-}) GeV", "MuMu_mass_X");

    compare_years("HLTemulation", "M_PiPi", "M_X3872 < 3.75", 50, 0., 1., "M(#pi^{+}#pi^{-}) GeV", "PiPi_mass_Psi");
    compare_years("HLTemulation", "M_PiPi", "M_X3872 > 3.75", 50, 0., 1., "M(#pi^{+}#pi^{-}) GeV", "PiPi_mass_X");

    compare_years("HLTemulation", "M_X3872", "M_X3872 < 3.75", 50, 3.50, 3.75,  "M(#psi (2S)) GeV", "Psi_mass_Psi");
    compare_years("HLTemulation", "M_X3872", "M_X3872 > 3.75", 50, 3.75, 4.0, "M(X_{3872}) GeV", "X_mass_X");

    compare_years("HLTemulation", "M_K0s", "M_X3872 < 3.75", 50, 0.4, 0.6,  "M(K_^{0}_{s}) GeV", "K0s_mass_Psi");
    compare_years("HLTemulation", "M_K0s", "M_X3872 > 3.75", 50, 0.4, 0.6, "M(K_^{0}_{s}) GeV", "K0s_mass_X");

    compare_years("HLTemulation", "pTM_B0", "", 15, 0., 15, "p_{T}(B^{0})/M_(B^{0})"); 
    compare_years("HLTemulation", "LxySignBSz_B0", "", 50, 0., 100, "L_{xy}/#sigma B^{0}-BS_{z}");
    compare_years("HLTemulation", "SVprob_B0", "", 10, 0., 1., "P_{sv}(B^{0})");
    gPad->SetLogy(1);
    compare_years("HLTemulation", "CosAlpha3DBSz_B0", "", 20, 0.999, 1., "cos #alpha_{3D} (B^{0}-BS_{z})");
    gPad->SetLogy(0);
    compare_years("HLTemulation", "LxySignSV_K0s", "", 50, 0., 100, "L_{xy}/#sigma B^{0}-K^{0}_{s}");
    compare_years("HLTemulation", "SVprob_PiPi", "", 10, 0., 1., "P_{sv}(#pi #pi)");
    compare_years("HLTemulation", "pT_PiPi", "", 40, 0., 0.4, "p_{T}(#pi #pi)/p_{T}(B^{0})");
    compare_years("HLTemulation", "pT_Pi1", "", 30, 0., 0.3, "p_{T}(#pi_{1})/p_{T}(B^{0})");
    compare_years("HLTemulation", "D0_Pi1", "", 30, 0., 30, "DCA(#pi_{1}, BS)");
    compare_years("HLTemulation", "DR_B0Pi1", "", 60, 0., 0.6, "#Delta R(#pi_{1}, B^{0})");
   


}
