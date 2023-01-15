#include "../include/RecoDecayX.h"

RecoDecayX::RecoDecayX(TTree *tree, const TString & tags) : MCbase_B0toX3872K0s (tree, tags){

    RecoP4_Mu1.SetM(mMuon); RecoP4_Mu2.SetM(mMuon);
    RecoP4_Pi1.SetM(mPion); RecoP4_Pi2.SetM(mPion);
    RecoP4_K0s.SetM(mK0s);

    outFilePath_ = "./outRoot/RecoDecay_X3872_" + tags_ + ".root";


}//RecoDecayX()

RecoDecayX::~RecoDecayX(){

    delete outFile_;
};

void RecoDecayX::Loop(){

    Long64_t nentries = fChain->GetEntriesFast();
    const Long64_t Nbreak = nentries + 10;
    const Long64_t Nprint = (int)(nentries/20.);


    // **** HISTOGRAMS **** //
    TH1F h_Mu_SoftID_MC = TH1F("Mu_SoftID_MC", "", 2, -0.5, 1.5);
    TH1F h_Mu_SoftID_Fk = TH1F("Mu_SoftID_Fk", "", 2, -0.5, 1.5);
    TH1F h_Mu_GlobalMu_MC = TH1F("Mu_GlobalMu_MC", "", 2, -0.5, 1.5);
    TH1F h_Mu_GlobalMu_Fk = TH1F("Mu_GlobalMu_Fk", "", 2, -0.5, 1.5);
    TH1F h_Mu_TrkQlty_MC = TH1F("Mu_TrkQlty_MC", "", 3, -0.5, 2.5);
    TH1F h_Mu_TrkQlty_Fk = TH1F("Mu_TrkQlty_Fk", "", 3, -0.5, 2.5);

    TH1F h_Pi_TrkQlty_MC = TH1F("Pi_TrkQlty_MC", "", 3, -0.5, 2.5);
    TH1F h_Pi_TrkQlty_Fk = TH1F("Pi_TrkQlty_Fk", "", 3, -0.5, 2.5);

    int Nbins = 20;
    double xlow = 0., xhigh = 0.005;
    TH1F h_Mu_dR_HLT_Dimuon25_Jpsi_MC = TH1F("Mu_dR_HLT_Dimuon25_Jpsi_MC", "", Nbins, xlow, xhigh);
    TH1F h_Mu_dR_HLT_Dimuon25_Jpsi_Fk = TH1F("Mu_dR_HLT_Dimuon25_Jpsi_Fk", "", Nbins, xlow, xhigh);
    TH1F h_Mu_dR_HLT_DoubleMu4_JpsiTrk_MC = TH1F("Mu_dR_HLT_DoubleMu4_JpsiTrk_MC", "", Nbins, xlow, xhigh);
    TH1F h_Mu_dR_HLT_DoubleMu4_JpsiTrk_Fk = TH1F("Mu_dR_HLT_DoubleMu4_JpsiTrk_Fk", "", Nbins, xlow, xhigh);
    xlow = 0., xhigh = 0.1;
    TH1F h_Pi_dR_HLT_DoubleMu4_JpsiTrk_MC = TH1F("Pi_dR_HLT_DoubleMu4_JpsiTrk_MC", "", Nbins, xlow, xhigh);
    TH1F h_Pi_dR_HLT_DoubleMu4_JpsiTrk_Fk = TH1F("Pi_dR_HLT_DoubleMu4_JpsiTrk_Fk", "", Nbins, xlow, xhigh);
    Nbins = 50, xlow = 0., xhigh = 1.;
    TH1F h_PiPi_svProb_MC = TH1F("PiPi_svProb_MC", "", Nbins, xlow, xhigh);
    TH1F h_PiPi_svProb_Fk = TH1F("PiPi_svProb_Fk", "", Nbins, xlow, xhigh);

    xlow = 0., xhigh = 100;
    TH1F h_K0s_LxySign_wrtBvtx_MC = TH1F("K0s_LxySign_wrtBvtx_MC", "", Nbins, xlow, xhigh);
    TH1F h_K0s_LxySign_wrtBvtx_Fk = TH1F("K0s_LxySign_wrtBvtx_Fk", "", Nbins, xlow, xhigh);
    Nbins = 20, xlow = 0., xhigh = 1.;
    TH1F h_K0s_cosAlpha2D_MC = TH1F("K0s_cosAlpha2D_MC", "", Nbins, xlow, xhigh);
    TH1F h_K0s_cosAlpha2D_Fk = TH1F("K0s_cosAlpha2D_Fk", "", Nbins, xlow, xhigh);
    TH1F h_K0s_cosAlpha3D_MC = TH1F("K0s_cosAlpha3D_MC", "", Nbins, xlow, xhigh);
    TH1F h_K0s_cosAlpha3D_Fk = TH1F("K0s_cosAlpha3D_Fk", "", Nbins, xlow, xhigh);
    xlow = 0.999, xhigh = 1.;
    TH1F h_B0_cosAlpha2D_MC = TH1F("B0_cosAlpha2D_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_cosAlpha2D_Fk = TH1F("B0_cosAlpha2D_Fk", "", Nbins, xlow, xhigh);
    TH1F h_B0_cosAlpha2DwrtBS_MC = TH1F("B0_cosAlpha2DwrtBS_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_cosAlpha2DwrtBS_Fk = TH1F("B0_cosAlpha2DwrtBS_Fk", "", Nbins, xlow, xhigh);
    TH1F h_B0_cosAlpha3D_MC = TH1F("B0_cosAlpha3D_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_cosAlpha3D_Fk = TH1F("B0_cosAlpha3D_Fk", "", Nbins, xlow, xhigh);
    Nbins = 50, xlow = 0., xhigh = 100;
    TH1F h_B0_LxySign_wrtPV_MC = TH1F("B0_LxySign_wrtPV_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_LxySign_wrtPV_Fk = TH1F("B0_LxySign_wrtPV_Fk", "", Nbins, xlow, xhigh);
    TH1F h_B0_LxySign_wrtBS_MC = TH1F("B0_LxySign_wrtBS_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_LxySign_wrtBS_Fk = TH1F("B0_LxySign_wrtBS_Fk", "", Nbins, xlow, xhigh);


    // MCmatching variables
    bool isMCmatched_Mu1, isMCmatched_Mu2, isMCmatched_Pi1, isMCmatched_Pi2;
    bool isMCmatched_JPsi, isMCmatched_Rho, isMCmatched_X3872, isMCmatched_K0s, isMCmatched_B0;
    int N_B0matching = 0;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0 || jentry == Nbreak) break;
        if ((jentry+1) % Nprint == 0) std::cout << "--> " << Form("%3.0f",(float)(jentry+1)/nentries* 100.) << " \%"<< std::endl;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        // ----- GENERATOR
        GenPartFillP4();

        // ----- FIND THE MONTE CARLO TRUTH
        //std::cout << " --- EV " << jentry << " with #B0 candidates " << nB0 << std::endl;
        MCtruthMatching();

        for (Int_t b = 0; b  < nB0; b++){
            
           if ( !RecoPartFillP4(b)) continue; 
            
        // --> check the MC matching of the chain
        isMCmatched_Mu1   =  (B0_mu1_idx[b] == MCmatch_Mum_Idx) || (B0_mu1_idx[b] == MCmatch_Mup_Idx);
        isMCmatched_Mu2   =  (B0_mu2_idx[b] == MCmatch_Mum_Idx) || (B0_mu2_idx[b] == MCmatch_Mup_Idx);
        isMCmatched_JPsi  =  isMCmatched_Mu1 && isMCmatched_Mu2;
        isMCmatched_Pi1   =  (B0_pi1_idx[b] == MCmatch_Pim_Idx) || (B0_pi1_idx[b] == MCmatch_Pip_Idx);
        isMCmatched_Pi2   =  (B0_pi2_idx[b] == MCmatch_Pim_Idx) || (B0_pi2_idx[b] == MCmatch_Pip_Idx);
        isMCmatched_Rho   =  isMCmatched_Pi1 && isMCmatched_Pi2;
        isMCmatched_X3872 =  isMCmatched_Rho && isMCmatched_JPsi;
        isMCmatched_K0s   =  (fabs(ROOT::Math::VectorUtil::DeltaR(GenP4_K0s, RecoP4_K0s) - MCmatch_K0s_DRmin) < 0.0001) && (MCmatch_K0s_Idx  != -1);
        isMCmatched_B0    =  isMCmatched_X3872 &&   isMCmatched_K0s;   



        // *** JPsi --> Mu+Mu-  ***
        if (isMCmatched_Mu1){
            h_Mu_SoftID_MC.Fill(Muon_softId[B0_mu1_idx[b]]);
            h_Mu_GlobalMu_MC.Fill(Muon_isGlobal[B0_mu1_idx[b]]);
            if (Muon_trackQuality[B0_mu1_idx[b]] < 2) h_Mu_TrkQlty_MC.Fill(Muon_trackQuality[B0_mu1_idx[b]]);
            else h_Mu_TrkQlty_MC.Fill(2);

            h_Mu_dR_HLT_Dimuon25_Jpsi_MC.Fill(B0_MuMu_mu1_dr_Dimuon25_Jpsi[b]);
            h_Mu_dR_HLT_DoubleMu4_JpsiTrk_MC.Fill(B0_MuMu_mu1_dr_DoubleMu4_JpsiTrk_Displaced[b]);
        }else{
            h_Mu_SoftID_Fk.Fill(Muon_softId[B0_mu1_idx[b]]);
            h_Mu_GlobalMu_Fk.Fill(Muon_isGlobal[B0_mu1_idx[b]]);
            if (Muon_trackQuality[B0_mu1_idx[b]] < 2) h_Mu_TrkQlty_Fk.Fill(Muon_trackQuality[B0_mu1_idx[b]]);
            else h_Mu_TrkQlty_Fk.Fill(2);

            h_Mu_dR_HLT_Dimuon25_Jpsi_Fk.Fill(B0_MuMu_mu1_dr_Dimuon25_Jpsi[b]);
            h_Mu_dR_HLT_DoubleMu4_JpsiTrk_Fk.Fill(B0_MuMu_mu1_dr_DoubleMu4_JpsiTrk_Displaced[b]);
        }
        if (isMCmatched_Mu2){
            h_Mu_SoftID_MC.Fill(Muon_softId[B0_mu2_idx[b]]);
            h_Mu_GlobalMu_MC.Fill(Muon_isGlobal[B0_mu2_idx[b]]);
            if (Muon_trackQuality[B0_mu1_idx[b]] < 2) h_Mu_TrkQlty_MC.Fill(Muon_trackQuality[B0_mu2_idx[b]]);
            else  h_Mu_TrkQlty_MC.Fill(2);

            h_Mu_dR_HLT_Dimuon25_Jpsi_MC.Fill(B0_MuMu_mu2_dr_Dimuon25_Jpsi[b]);
            h_Mu_dR_HLT_DoubleMu4_JpsiTrk_MC.Fill(B0_MuMu_mu2_dr_DoubleMu4_JpsiTrk_Displaced[b]);
        }else{
            h_Mu_SoftID_Fk.Fill(Muon_softId[B0_mu2_idx[b]]);
            h_Mu_GlobalMu_Fk.Fill(Muon_isGlobal[B0_mu2_idx[b]]);
            if (Muon_trackQuality[B0_mu1_idx[b]] < 2) h_Mu_TrkQlty_Fk.Fill(Muon_trackQuality[B0_mu2_idx[b]]);
            else h_Mu_TrkQlty_Fk.Fill(2);

            h_Mu_dR_HLT_Dimuon25_Jpsi_Fk.Fill(B0_MuMu_mu2_dr_Dimuon25_Jpsi[b]);
            h_Mu_dR_HLT_DoubleMu4_JpsiTrk_Fk.Fill(B0_MuMu_mu2_dr_DoubleMu4_JpsiTrk_Displaced[b]);
        }
        
        
        // *** Rho --> Pi+Pi-  ***
        if(isMCmatched_Pi1){
            h_Pi_dR_HLT_DoubleMu4_JpsiTrk_MC.Fill(B0_PiPi_p1_dr_DoubleMu4_JpsiTrk_Displaced[b]);
        }else{
            h_Pi_dR_HLT_DoubleMu4_JpsiTrk_Fk.Fill(B0_PiPi_p1_dr_DoubleMu4_JpsiTrk_Displaced[b]);
            //std::cout << B0_PiPi_p1_dr_DoubleMu4_JpsiTrk_Displaced[b]<< std::endl;
        }
        if (isMCmatched_Pi2)
        {
            h_Pi_dR_HLT_DoubleMu4_JpsiTrk_MC.Fill(B0_PiPi_p2_dr_DoubleMu4_JpsiTrk_Displaced[b]);   
        }else{
            h_Pi_dR_HLT_DoubleMu4_JpsiTrk_Fk.Fill(B0_PiPi_p2_dr_DoubleMu4_JpsiTrk_Displaced[b]);
        }
        if (isMCmatched_Rho){
            h_PiPi_svProb_MC.Fill(B0_PiPi_sv_prob[b]);
        }else{
            h_PiPi_svProb_Fk.Fill(B0_PiPi_sv_prob[b]);
        }
        

        // *** K0short ***
        if (isMCmatched_K0s){
            h_K0s_LxySign_wrtBvtx_MC.Fill(B0_K0_lxySign_wrtBvtx[b]);
            h_K0s_cosAlpha2D_MC.Fill(fabs(B0_K0_cosAlpha2D[b]));
            h_K0s_cosAlpha3D_MC.Fill(fabs(B0_K0_cosAlpha3D[b]));
 
        }else{
            h_K0s_LxySign_wrtBvtx_Fk.Fill(B0_K0_lxySign_wrtBvtx[b]);
            h_K0s_cosAlpha2D_Fk.Fill(fabs(B0_K0_cosAlpha2D[b]));
            h_K0s_cosAlpha3D_Fk.Fill(fabs(B0_K0_cosAlpha3D[b]));

        }
        
        
        // *** B0 --> X(3872) K0s ***
        if (isMCmatched_B0){
            h_B0_cosAlpha2D_MC.Fill(fabs(B0_cosAlpha2D_PV[b]));
            h_B0_cosAlpha2DwrtBS_MC.Fill(fabs(B0_cosAlpha2D_BS[b]));
            h_B0_cosAlpha3D_MC.Fill(fabs(B0_cosAlpha3D_PV[b]));

            h_B0_LxySign_wrtPV_MC.Fill(B0_lxySign_PV[b]);
            h_B0_LxySign_wrtBS_MC.Fill(B0_lxySign_BS[b]);
        }else{
            h_B0_cosAlpha2D_Fk.Fill(fabs(B0_cosAlpha2D_PV[b]));
            h_B0_cosAlpha2DwrtBS_Fk.Fill(fabs(B0_cosAlpha2D_BS[b]));
            h_B0_cosAlpha3D_Fk.Fill(fabs(B0_cosAlpha3D_PV[b]));

            h_B0_LxySign_wrtPV_Fk.Fill(B0_lxySign_PV[b]);
            h_B0_LxySign_wrtBS_Fk.Fill(B0_lxySign_BS[b]);
        }

        if (isMCmatched_B0) N_B0matching++;

        
        } // loop on B0 candidates
       


    }// loop on events

    std::cout << " B0 cand. MC matching " << N_B0matching << std::endl;
    
    outFile_ = new TFile(outFilePath_, "RECREATE");
    outFile_->cd();
    h_Mu_SoftID_MC.Write();
    h_Mu_SoftID_Fk.Write();
    h_Mu_GlobalMu_MC.Write();
    h_Mu_GlobalMu_Fk.Write();
    h_Mu_TrkQlty_MC.Write();
    h_Mu_TrkQlty_Fk.Write();

    h_Mu_dR_HLT_Dimuon25_Jpsi_MC.Write();
    h_Mu_dR_HLT_Dimuon25_Jpsi_Fk.Write();
    h_Mu_dR_HLT_DoubleMu4_JpsiTrk_MC.Write();
    h_Mu_dR_HLT_DoubleMu4_JpsiTrk_Fk.Write();

    h_Pi_dR_HLT_DoubleMu4_JpsiTrk_MC.Write();
    h_Pi_dR_HLT_DoubleMu4_JpsiTrk_Fk.Write();
    h_PiPi_svProb_MC.Write();
    h_PiPi_svProb_Fk.Write();

    h_K0s_LxySign_wrtBvtx_MC.Write();
    h_K0s_LxySign_wrtBvtx_Fk.Write();
    h_K0s_cosAlpha2D_MC.Write();
    h_K0s_cosAlpha2D_Fk.Write();
    h_K0s_cosAlpha3D_MC.Write();
    h_K0s_cosAlpha3D_Fk.Write();


    h_B0_cosAlpha2D_MC.Write();
    h_B0_cosAlpha2D_Fk.Write();
    h_B0_cosAlpha2DwrtBS_MC.Write();
    h_B0_cosAlpha2DwrtBS_Fk.Write();
    h_B0_cosAlpha3D_MC.Write();
    h_B0_cosAlpha3D_Fk.Write();
    h_B0_LxySign_wrtPV_MC.Write();
    h_B0_LxySign_wrtPV_Fk.Write();
    h_B0_LxySign_wrtBS_MC.Write();
    h_B0_LxySign_wrtBS_Fk.Write();

    outFile_->Close();
    std::cout << "  ...[OUTPUT] output histograms written on file " << outFilePath_ << std::endl;


}//Loop()


int RecoDecayX::RecoPartFillP4(const int Bidx){
    int TrackQualityCheck = 1;
 
    //... muons P4
    if(!Muon_softId[B0_mu1_idx[Bidx]] || !Muon_softId[B0_mu2_idx[Bidx]]) TrackQualityCheck = 0;
    RecoP4_Mu1.SetPt(B0_MuMu_prefit_mu1_pt[Bidx]); RecoP4_Mu1.SetEta(B0_MuMu_prefit_mu1_eta[Bidx]); RecoP4_Mu1.SetPhi(B0_MuMu_prefit_mu1_phi[Bidx]);
    RecoP4_Mu2.SetPt(B0_MuMu_prefit_mu2_pt[Bidx]); RecoP4_Mu2.SetEta(B0_MuMu_prefit_mu2_eta[Bidx]); RecoP4_Mu2.SetPhi(B0_MuMu_prefit_mu2_phi[Bidx]);
    //... tracks P4
    if(ProbeTracks_isMatchedToMuon[B0_pi1_idx[Bidx]] || ProbeTracks_isMatchedToMuon[B0_pi2_idx[Bidx]])TrackQualityCheck = 0;
    RecoP4_Pi1.SetPt(B0_PiPi_prefit_pi1_pt[Bidx]); RecoP4_Pi1.SetEta(B0_PiPi_prefit_pi1_eta[Bidx]); RecoP4_Pi1.SetPhi(B0_PiPi_prefit_pi1_phi[Bidx]);
    RecoP4_Pi2.SetPt(B0_PiPi_prefit_pi2_pt[Bidx]); RecoP4_Pi2.SetEta(B0_PiPi_prefit_pi2_eta[Bidx]); RecoP4_Pi2.SetPhi(B0_PiPi_prefit_pi2_phi[Bidx]);
    RecoP4_K0s.SetPt(B0_K0s_mcFitted_pt[Bidx]); RecoP4_K0s.SetEta(B0_K0s_mcFitted_eta[Bidx]); RecoP4_K0s.SetPhi(B0_K0s_mcFitted_phi[Bidx]);
 
    return TrackQualityCheck;
 
 }//RecoPartFillP4()

void RecoDecayX::MCtruthMatching(const bool verbose){

    const double DRmin_threshold = 0.03;
    const double DpT_threshold = 0.5;
    float DeltaPt;

    // ..... muons ..... //	
    float DRminMum = 100., DR_gMum_rMu1, DR_gMum_rMu2, DRminMum_DpT = 0.; 
    float DRminMup = 100., DR_gMup_rMu1, DR_gMup_rMu2, DRminMup_DpT = 0.;
    MCmatch_Mum_Idx = -1, MCmatch_Mup_Idx = -1; 
    // ..... pions ..... //      
    float DRminPim = 100., DR_gPim_rPi1, DR_gPim_rPi2, DRminPim_DpT = 0.;
    float DRminPip = 100., DR_gPip_rPi1, DR_gPip_rPi2, DRminPip_DpT = 0.;
    MCmatch_Pim_Idx = -1, MCmatch_Pip_Idx = -1;

    // ..... K0s ..... //
    float DRminK0s = 100., DR_gK0s_rK0s, DRminK0s_DpT =0;
    MCmatch_K0s_Idx = -1;

   for (UInt_t b = 0; b < nB0; b++){

       RecoPartFillP4(b);

      // ..... muons ..... //
		if (Muon_softId[B0_mu1_idx[b]] && Muon_softId[B0_mu2_idx[b]]){ // softId for both muons is required

			DR_gMum_rMu1 = ROOT::Math::VectorUtil::DeltaR(GenP4_Mum, RecoP4_Mu1); // mu1
			DR_gMup_rMu1 = ROOT::Math::VectorUtil::DeltaR(GenP4_Mup, RecoP4_Mu1);
			DR_gMum_rMu2 = ROOT::Math::VectorUtil::DeltaR(GenP4_Mum, RecoP4_Mu2); // mu2
			DR_gMup_rMu2 = ROOT::Math::VectorUtil::DeltaR(GenP4_Mup, RecoP4_Mu2);

			// DR(mu-, mu1) < DRmin(mu-) + m1 is nearer to mu- than mu2 + DRmin threshold
			// mu-
			if( (DR_gMum_rMu1 < DRminMum) && (DR_gMum_rMu1 < DR_gMum_rMu2) ){
				DRminMum = DR_gMum_rMu1;
				DRminMum_DpT = DeltaPT(GenP4_Mum, RecoP4_Mu1);
				MCmatch_Mum_Idx = B0_mu1_idx[b];
			}
			if( (DR_gMum_rMu2 < DRminMum) && (DR_gMum_rMu2 < DR_gMum_rMu1) ){
				DRminMum = DR_gMum_rMu2;
				DRminMum_DpT = DeltaPT(GenP4_Mum, RecoP4_Mu2);
				MCmatch_Mum_Idx = B0_mu2_idx[b];
			}
			// mu+
			if( (DR_gMup_rMu1 < DRminMup) && (DR_gMup_rMu1 < DR_gMup_rMu2) ){
				DRminMup = DR_gMup_rMu1;
				DRminMup_DpT = DeltaPT(GenP4_Mup, RecoP4_Mu1);
				MCmatch_Mup_Idx = B0_mu1_idx[b];
			}
			if( (DR_gMup_rMu2 < DRminMup) && (DR_gMup_rMu2 < DR_gMup_rMu1) ){
				DRminMup = DR_gMup_rMu2;
				DRminMup_DpT = DeltaPT(GenP4_Mup, RecoP4_Mu2);
				MCmatch_Mup_Idx = B0_mu2_idx[b];
			}

		}// ... muons //

		// ..... pions ..... //
		if (!ProbeTracks_isMatchedToMuon[B0_pi1_idx[b]] && !ProbeTracks_isMatchedToMuon[B0_pi2_idx[b]]){ // isMatchedToMuon must be false

			DR_gPim_rPi1 = ROOT::Math::VectorUtil::DeltaR(GenP4_Pim, RecoP4_Pi1);//pi1
			DR_gPip_rPi1 = ROOT::Math::VectorUtil::DeltaR(GenP4_Pip, RecoP4_Pi1);
			DR_gPim_rPi2 = ROOT::Math::VectorUtil::DeltaR(GenP4_Pim, RecoP4_Pi2);//pi2
			DR_gPip_rPi2 = ROOT::Math::VectorUtil::DeltaR(GenP4_Pip, RecoP4_Pi2);
			// pi-
			if( (DR_gPim_rPi1 < DRminPim) && (DR_gPim_rPi1 < DR_gPim_rPi2) ){
				DRminPim = DR_gPim_rPi1;
				DRminPim_DpT = DeltaPT(GenP4_Pim, RecoP4_Pi1);
				MCmatch_Pim_Idx = B0_pi1_idx[b];
			}
			if( (DR_gPim_rPi2 < DRminPim) && (DR_gPim_rPi2 < DR_gPim_rPi1) ){
				DRminPim = DR_gPim_rPi2;
				DRminPim_DpT = DeltaPT(GenP4_Pim, RecoP4_Pi2);
				MCmatch_Pim_Idx = B0_pi2_idx[b];
			}
			// pi+
			if( (DR_gPip_rPi1 < DRminPip) && (DR_gPip_rPi1 < DR_gPip_rPi2) ){
				DRminPip = DR_gPip_rPi1;
				DRminPip_DpT = DeltaPT(GenP4_Pip, RecoP4_Pi1);
				MCmatch_Pip_Idx = B0_pi1_idx[b];
			}
			if( (DR_gPip_rPi2 < DRminPip) && (DR_gPip_rPi2 < DR_gPip_rPi1) ){
				DRminPip = DR_gPip_rPi2;
				DRminPip_DpT = DeltaPT(GenP4_Pip, RecoP4_Pi2);
				MCmatch_Pip_Idx = B0_pi2_idx[b];
			}


		}// ... pions//

		// ..... K0s ..... //
		DR_gK0s_rK0s = ROOT::Math::VectorUtil::DeltaR(GenP4_K0s, RecoP4_K0s);
		if ( (DR_gK0s_rK0s < DRminK0s) ){ 
			DRminK0s = DR_gK0s_rK0s;
			DRminK0s_DpT = DeltaPT(GenP4_K0s, RecoP4_K0s);
			MCmatch_K0s_Idx = B0_k0short_idx[b];
		}


	}// on B0 candidates

    // ... muons
	if(DRminMum > DRmin_threshold) MCmatch_Mum_Idx = -1; 
	MCmatch_Mum_DRmin = DRminMum; MCmatch_Mum_DpT = DRminMum_DpT; 
	if(DRminMup > DRmin_threshold) MCmatch_Mup_Idx = -1; 
	MCmatch_Mup_DRmin = DRminMup; MCmatch_Mup_DpT = DRminMup_DpT; 

    // ... pions
	if( (DRminPim > DRmin_threshold) || (DRminPim_DpT > DpT_threshold)) MCmatch_Pim_Idx = -1;
	MCmatch_Pim_DRmin = DRminPim; MCmatch_Pim_DpT = DRminPim_DpT;
	if( (DRminPip > DRmin_threshold) || (DRminPip_DpT > DpT_threshold)) MCmatch_Pip_Idx = -1;
	MCmatch_Pip_DRmin = DRminPip; MCmatch_Pip_DpT = DRminPip_DpT; 

    // ... K0 short
	if( (DRminK0s > DRmin_threshold) || (DRminK0s_DpT > DpT_threshold)) MCmatch_K0s_Idx = -1;
	MCmatch_K0s_DRmin = DRminK0s; MCmatch_K0s_DpT = DRminK0s_DpT;
 
	if (verbose){
		std::cout << "MC matching indices " << std::endl;
		std::cout << "mu- " << MCmatch_Mum_Idx << "\t mu+ " << MCmatch_Mup_Idx << std::endl;
        //std::cout << "DR  " << MCmatch_Mum_DRmin << "\t mu+ " << MCmatch_Mup_DRmin << std::endl;
		std::cout << "pi- " << MCmatch_Pim_Idx << "\t pi+ " << MCmatch_Pip_Idx << std::endl;
        //std::cout << "DR  " << MCmatch_Pim_DRmin << "\t mu+ " << MCmatch_Pip_DRmin << std::endl;
		std::cout << "K0s " << MCmatch_K0s_Idx << std::endl;
        //std::cout << "DR  " << MCmatch_K0s_DRmin << std::endl;
	}

}//MCtruthMatching()

float RecoDecayX::DeltaPT(ROOT::Math::PtEtaPhiMVector& genV, ROOT::Math::PtEtaPhiMVector& recV){
	return fabs(genV.Pt() - recV.Pt()) / genV.Pt();
}//DeltaPT()
