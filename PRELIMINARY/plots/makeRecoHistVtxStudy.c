{
    //gStyle->SetLegendTextSize(0.03);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    //.L plots/histoManagerReco.cc++;

    SetInputFile("./outRoot/RecoDecay_X3872_UL17_TrkFixVtxStudy.root");
    SetOutputFile("/eos/user/c/cbasile/www/B0toX3872K0s/RECO_LEVEL/RecoDecay_X3872_UL17_TrkFixVtxStudy/");

    draw_one_histo("Mu_dR_HLT_Dimuon25_Jpsi_MC", "mu_MC", "#Delta R (#mu_{HLT}, #mu_{1/2})", "");
    draw_one_histo("Mu_dR_HLT_DoubleMu4_JpsiTrk_MC", "mu_MC", "#Delta R (#mu_{HLT}, #mu_{1/2})", "");
    draw_binary_histo("Mu_SoftID_MC", "mu_MC", "muon SoftID","Mu_SoftID");
    draw_binary_histo("Mu_GlobalMu_MC", "mu_MC", "Global muon","Mu_isGlobal");
    draw_binary_histo("Mu_TrkQlty_MC", "mu_MC", "muon Tarck Quality","Mu_TrkQuality");

    draw_two_histograms("Pi_dR_HLT_DoubleMu4_JpsiTrk_MC", "pi_MC", "Pi_dR_HLT_DoubleMu4_JpsiTrk_Fk", "pi_Fk", 
            "#Delta R (trk_{HLT}, trk_{reco})", "Pi_dR_HLT_DoubleMu4_JpsiTrk");
    draw_two_histograms("PiPi_svProb_MC", "MCmatch", "PiPi_svProb_Fk", "MCfake", "P_{sv}(#pi_{1} #pi_{2})", "PiPi_svProb" );
     
    
    draw_two_histograms("K0s_LxySign_wrtBvtx_MC", "MCmatch", "K0s_LxySign_wrtBvtx_Fk", "MCfake", 
            "L_{xy}(K^{0_s}; B^{0}_{vtx})/#sigma_{xy}", "K0s_LxySign_wrtBvtx" );

    draw_two_histograms("B0_cosAlpha3DwrtPV_MC", "MCmatch", "B0_cosAlpha3DwrtPV_Fk", "MCfake", 
            "cos #alpha (#vec{p}_{B}, L_{B})_{3D}", "B0_cosAlpha3DwrtPV", true, true);
    draw_two_histograms("B0_cosAlpha2DwrtBSrk_MC", "MCmatch", "B0_cosAlpha2DwrtBSrk_Fk", "MCfake", 
            "cos #alpha (#vec{p}_{B}, L_{B})_{2D}", "B0_cosAlpha2DwrtBSrk", true, true);
    draw_two_histograms("B0_cosAlpha2DwrtBS_MC", "MCmatch", "B0_cosAlpha2DwrtBS_Fk", "MCfake", 
            "cos #alpha (#vec{p}_{B}, L_{B})_{2D}", "B0_cosAlpha2DwrtBS", true, true);
    
    
    draw_two_histograms("B0_LxySign_wrtPV_MC", "MCmatch", "B0_LxySign_wrtPV_Fk", "MCfake", 
            "L_{xy}( PV, B^{0}_{vtx})/#sigma_{xy}", "B0_LxySign_wrtPV" );
    draw_two_histograms("B0_LxySign_wrtBS_MC", "MCmatch", "B0_LxySign_wrtBS_Fk", "MCfake", 
            "L_{xy}( BS, B^{0}_{vtx})/#sigma_{xy}", "B0_LxySign_wrtBS" );
    draw_two_histograms("B0_LxySign_wrtBSrk_MC", "MCmatch", "B0_LxySign_wrtBSrk_Fk", "MCfake", 
            "L_{xy}( BS, B^{0}_{vtx})/#sigma_{xy}", "B0_LxySign_wrtBSrk" );
    draw_two_histograms("B0_LxySign_wrtBShlt_MC", "MCmatch", "B0_LxySign_wrtBShlt_Fk", "MCfake", 
            "L_{xy}( BS, B^{0}_{vtx})/#sigma_{xy}", "B0_LxySign_wrtBShlt" );
    
    draw_two_histograms("B0_DistGenVtx_PV_MC", "PV", "B0_DistGenVtx_BS_MC", "BS", "|V_{G} - V_{R}|", "DistGenVTXwrtReco" );


// ROC curves
    makeROCcurve({"B0_cosAlpha2DwrtBS_MC","B0_cosAlpha3DwrtPV_MC","B0_cosAlpha2DwrtBSrk_MC"},{"B0_cosAlpha2DwrtBS_Fk","B0_cosAlpha3DwrtPV_Fk","B0_cosAlpha2DwrtBSrk_Fk"}, "ROC_B0_cosAlpha");

    makeROCcurve({"B0_LxySign_wrtPV_MC", "B0_LxySign_wrtBS_MC"},{"B0_LxySign_wrtPV_Fk", "B0_LxySign_wrtBS_Fk"}, "ROC_LxySign_B0_PVvsBS");
    makeROCcurve({"B0_LxySign_wrtBS_MC", "B0_LxySign_wrtBShlt_MC","B0_LxySign_wrtBSrk_MC"}, {"B0_LxySign_wrtBS_Fk", "B0_LxySign_wrtBShlt_Fk","B0_LxySign_wrtBSrk_Fk"}, "ROC_LxySign_B0wrtBS");
    makeROCcurve({"B0_LxySign_wrtPV_MC", "B0_LxySign_wrtBSrk_MC"},{"B0_LxySign_wrtPV_Fk", "B0_LxySign_wrtBSrk_Fk"}, "ROC_LxySign_B0_PVvsBSrk");

}
