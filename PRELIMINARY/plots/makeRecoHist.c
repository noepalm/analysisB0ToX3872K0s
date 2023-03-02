{
    //gStyle->SetLegendTextSize(0.03);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    //.L plots/histoManagerReco.cc++;

    SetInputFile("./outRoot/RecoDecay_X3872_UL16preVFP.root");
    SetOutputFile("/eos/user/c/cbasile/www/B0toX3872K0s/RECO_LEVEL/RecoDecay_X3872_UL16preVFP/");


    // MASS PLOTS
    draw_one_histo("MuMu_M_MC", "mu_MC", "M(mu_{1} #mu_{2})", "MuMu_mass");
    draw_two_histograms("PiPi_M_MC", "PiPi", "PiPi_M_Fk", "PiPi_Fk", "M(#pi_{1} #pi_{2})", "PiPi_mass" );
    draw_two_histograms("X3872_M_MC", "X", "X3872_M_Fk", "X_Fk", "M(X(3872))", "X3872_mass" );
    draw_two_histograms("K0s_M_MC", "K0s", "K0s_M_Fk", "K0s_Fk", "M(K^{0}_s)", "K0s_mass" );
    draw_two_histograms("B0_M_MC", "MCmatch", "B0_M_Fk", "MCfake", "M(B^{0})", "B0_mass" );



    draw_one_histo("Mu_dR_HLT_Dimuon25_Jpsi_MC", "mu_MC", "#Delta R (#mu_{HLT}, #mu_{1/2})", "");
    draw_one_histo("Mu_dR_HLT_DoubleMu4_JpsiTrk_MC", "mu_MC", "#Delta R (#mu_{HLT}, #mu_{1/2})", "");
    draw_binary_histo("Mu_SoftID_MC", "mu_MC", "muon SoftID","Mu_SoftID");
    draw_binary_histo("Mu_GlobalMu_MC", "mu_MC", "Global muon","Mu_isGlobal");
    draw_binary_histo("Mu_TrkQlty_MC", "mu_MC", "muon Tarck Quality","Mu_TrkQuality");

    draw_two_histograms("Pi_dR_HLT_DoubleMu4_JpsiTrk_MC", "pi_MC", "Pi_dR_HLT_DoubleMu4_JpsiTrk_Fk", "pi_Fk", 
            "#Delta R (trk_{HLT}, trk_{reco})", "Pi_dR_HLT_DoubleMu4_JpsiTrk");
    
    draw_two_histograms("PiPi_svProb_MC", "MCmatch", "PiPi_svProb_Fk", "MCfake", "P_{sv}(#pi_{1} #pi_{2})", "PiPi_svProb" );
    draw_two_histograms("PiPi_pT_MC", "MCmatch", "PiPi_pT_Fk", "MCfake", "p_{T}(#pi #pi)/p_{T}(B^{0})", "PiPi_pT" );
    draw_two_histograms("Pi1_pT_MC", "MCmatch", "Pi1_pT_Fk", "MCfake", "p_{T}(#pi_{1})/p_{T}(B^{0})", "Pi1_pT" );
    draw_two_histograms("Pi1_D0_MC", "MCmatch", "Pi1_D0_Fk", "MCfake", "DCA(#pi_{1}, BS)", "Pi1_D0" );
    draw_two_histograms("Pi1_DRwrtB_MC", "MCmatch", "Pi1_DRwrtB_Fk", "MCfake", "#Delta R(#pi_{1}, B^{0})", "Pi1_DRwrtB0" );
     
    
    draw_two_histograms("K0s_LxySign_wrtBvtx_MC", "MCmatch", "K0s_LxySign_wrtBvtx_Fk", "MCfake", 
            "L_{xy}(K^{0_s}; B^{0}_{vtx})/#sigma_{xy}", "K0s_LxySign_wrtBvtx" );

    draw_two_histograms("B0_cosAlpha3DwrtPV_MC", "MCmatch", "B0_cosAlpha3DwrtPV_Fk", "MCfake", 
            "cos #alpha (#vec{p}_{B}, L_{B})_{3D}", "B0_cosAlpha3DwrtPV", true, true);
    draw_two_histograms("B0_cosAlpha2DwrtBSwithZ_MC", "MCmatch", "B0_cosAlpha2DwrtBSwithZ_Fk", "MCfake", 
            "cos #alpha (#vec{p}_{B}, L_{B})_{2D}", "B0_cosAlpha2DwrtBSwithZ", true, true);
    draw_two_histograms("B0_cosAlpha2DwrtBS_MC", "MCmatch", "B0_cosAlpha2DwrtBS_Fk", "MCfake", 
            "cos #alpha (#vec{p}_{B}, L_{B})_{2D}", "B0_cosAlpha2DwrtBS", true, true);
    
    
    draw_two_histograms("B0_LxySign_wrtPV_MC", "MCmatch", "B0_LxySign_wrtPV_Fk", "MCfake", 
            "L_{xy}( PV, B^{0}_{vtx})/#sigma_{xy}", "B0_LxySign_wrtPV" );
    draw_two_histograms("B0_LxySign_wrtBS_MC", "MCmatch", "B0_LxySign_wrtBS_Fk", "MCfake", 
            "L_{xy}( BS, B^{0}_{vtx})/#sigma_{xy}", "B0_LxySign_wrtBS" );
    draw_two_histograms("B0_LxySign_wrtBSwithZ_MC", "MCmatch", "B0_LxySign_wrtBSwithZ_Fk", "MCfake", 
            "L_{xy}( BS_{z}, B^{0}_{vtx})/#sigma_{xy}", "B0_LxySign_wrtBSwithZ" );
    
    draw_two_histograms("B0_svProb_MC", "MCmatch", "B0_svProb_Fk", "MCfake", "P_{sv}(B^{0})", "B0_svProb" );
    draw_two_histograms("B0_pT_MC", "MCmatch", "B0_pT_Fk", "MCfake", "p_{T}(B^{0})/M(B^{0})", "B0_pT" );
    
    draw_two_histograms("B0_DistGenVtx_PV_MC", "PV", "B0_DistGenVtx_BS_MC", "BS", "|V_{G} - V_{R}|", "DistGenVTXwrtReco" );


// ROC curves
    makeROCcurve({"B0_cosAlpha2DwrtBS_MC","B0_cosAlpha3DwrtPV_MC","B0_cosAlpha2DwrtBSwithZ_MC"},{"B0_cosAlpha2DwrtBS_Fk","B0_cosAlpha3DwrtPV_Fk","B0_cosAlpha2DwrtBSwithZ_Fk"}, "ROC_B0_cosAlpha");

    makeROCcurve({"B0_LxySign_wrtPV_MC", "B0_LxySign_wrtBS_MC", "B0_LxySign_wrtBSwithZ_MC"},{"B0_LxySign_wrtPV_Fk", "B0_LxySign_wrtBS_Fk", "B0_LxySign_wrtBSwithZ_Fk"}, "ROC_LxySign_B0_PVvsBS");
    makeROCcurve({"B0_LxySign_wrtPV_MC", "B0_LxySign_wrtBSwithZ_MC"},{"B0_LxySign_wrtPV_Fk", "B0_LxySign_wrtBSwithZ_Fk"}, "ROC_LxySign_B0_PVvsBSwithZ");

}
