#include "../include/RecoDecayPsi.h"

RecoDecayPsi::RecoDecayPsi(TTree *tree, const TString & tags) : MCbase_B0toPsi2SK0s (tree, tags){

    RecoP4_Mu1.SetM(mMuon); RecoP4_Mu2.SetM(mMuon);
    RecoP4_Pi1.SetM(mPion); RecoP4_Pi2.SetM(mPion);
    RecoP4_K0s.SetM(mK0s);

    outFilePath_ = "./outRoot/RecoDecay_Psi2S_" + tags_ + ".root";
    OutTree_setup();


}//RecoDecayPsi()

RecoDecayPsi::~RecoDecayPsi(){

    delete outFile_;
    delete outTree_;
}

void RecoDecayPsi::Loop(){

    Long64_t nentries = fChain->GetEntriesFast();
    const Long64_t Nbreak = nentries + 10; 
    const Long64_t Nprint = (int)(nentries/20.);

    int realB_idx;
    double Dist_GenV_PV, Dist_GenV_BS, Dist_GenV_BSwithZ;


    // **** HISTOGRAMS **** //
    // muons tracks quality
    TH1F h_Mu_SoftID_MC   = TH1F("Mu_SoftID_MC", "", 2, -0.5, 1.5);
    TH1F h_Mu_SoftID_Fk   = TH1F("Mu_SoftID_Fk", "", 2, -0.5, 1.5);
    TH1F h_Mu_GlobalMu_MC = TH1F("Mu_GlobalMu_MC", "", 2, -0.5, 1.5);
    TH1F h_Mu_GlobalMu_Fk = TH1F("Mu_GlobalMu_Fk", "", 2, -0.5, 1.5);
    TH1F h_Mu_TrkQlty_MC  = TH1F("Mu_TrkQlty_MC", "", 3, -0.5, 2.5);
    TH1F h_Mu_TrkQlty_Fk  = TH1F("Mu_TrkQlty_Fk", "", 3, -0.5, 2.5);
    // pions track quality
    TH1F h_Pi_TrkQlty_MC = TH1F("Pi_TrkQlty_MC", "", 3, -0.5, 2.5);
    TH1F h_Pi_TrkQlty_Fk = TH1F("Pi_TrkQlty_Fk", "", 3, -0.5, 2.5);

    int Nbins = 20;
    double xlow = 0., xhigh = 0.005;
    // DeltaR(trigger-Mu; reco-Mu)
    TH1F h_Mu_dR_HLT_Dimuon25_Jpsi_MC = TH1F("Mu_dR_HLT_Dimuon25_Jpsi_MC", "", Nbins, xlow, xhigh);
    TH1F h_Mu_dR_HLT_Dimuon25_Jpsi_Fk = TH1F("Mu_dR_HLT_Dimuon25_Jpsi_Fk", "", Nbins, xlow, xhigh);
    TH1F h_Mu_dR_HLT_DoubleMu4_JpsiTrk_MC = TH1F("Mu_dR_HLT_DoubleMu4_JpsiTrk_MC", "", Nbins, xlow, xhigh);
    TH1F h_Mu_dR_HLT_DoubleMu4_JpsiTrk_Fk = TH1F("Mu_dR_HLT_DoubleMu4_JpsiTrk_Fk", "", Nbins, xlow, xhigh);
    xlow = 0., xhigh = 0.03;
    // DeltaR(trigger-Pi; reco-Pi)
    TH1F h_Pi_dR_HLT_DoubleMu4_JpsiTrk_MC = TH1F("Pi_dR_HLT_DoubleMu4_JpsiTrk_MC", "", Nbins, xlow, xhigh);
    TH1F h_Pi_dR_HLT_DoubleMu4_JpsiTrk_Fk = TH1F("Pi_dR_HLT_DoubleMu4_JpsiTrk_Fk", "", Nbins, xlow, xhigh);
    
    // JPsi --> MuMu
    Nbins =35 , xlow = 2.9, xhigh = 3.25;
    TH1F h_MuMu_M_MC = TH1F("MuMu_M_MC", "", Nbins, xlow, xhigh);
    TH1F h_MuMu_M_Fk = TH1F("MuMu_M_Fk", "", Nbins, xlow, xhigh);

    // Rho --> PiPi
    Nbins = 30, xlow = 0., xhigh = .30; 
    TH1F h_Pi1_pT_MC        = TH1F("Pi1_pT_MC", "", Nbins, xlow, xhigh);
    TH1F h_Pi1_pT_Fk        = TH1F("Pi1_pT_Fk", "", Nbins, xlow, xhigh);
    Nbins = 30, xlow = 0., xhigh = 0.60; 
    TH1F h_Pi1_DRwrtB_MC    = TH1F("Pi1_DRwrtB_MC", "", Nbins, xlow, xhigh);
    TH1F h_Pi1_DRwrtB_Fk    = TH1F("Pi1_DRwrtB_Fk", "", Nbins, xlow, xhigh);
    Nbins = 15, xlow = 0., xhigh = 30; 
    TH1F h_Pi1_D0_MC        = TH1F("Pi1_D0_MC", "", Nbins, xlow, xhigh);
    TH1F h_Pi1_D0_Fk        = TH1F("Pi1_D0_Fk", "", Nbins, xlow, xhigh);
    Nbins = 50 , xlow = 0., xhigh = 1.;
    TH1F h_PiPi_M_MC        = TH1F("PiPi_M_MC", "", Nbins, xlow, xhigh);
    TH1F h_PiPi_M_Fk        = TH1F("PiPi_M_Fk", "", Nbins, xlow, xhigh);
    Nbins = 20, xlow = 0., xhigh = 1.;
    TH1F h_PiPi_svProb_MC   = TH1F("PiPi_svProb_MC", "", Nbins, xlow, xhigh);
    TH1F h_PiPi_svProb_Fk   = TH1F("PiPi_svProb_Fk", "", Nbins, xlow, xhigh);
    Nbins = 40, xlow = 0., xhigh = 0.4;
    TH1F h_PiPi_pT_MC       = TH1F("PiPi_pT_MC", "", Nbins, xlow, xhigh);
    TH1F h_PiPi_pT_Fk       = TH1F("PiPi_pT_Fk", "", Nbins, xlow, xhigh);

    // K0s --> PiPi
    xlow = 0., xhigh = 100;
    TH1F h_K0s_LxySign_wrtBvtx_MC = TH1F("K0s_LxySign_wrtBvtx_MC", "", Nbins, xlow, xhigh);
    TH1F h_K0s_LxySign_wrtBvtx_Fk = TH1F("K0s_LxySign_wrtBvtx_Fk", "", Nbins, xlow, xhigh);
    // B0 --> Psi K0s
    xlow = 0.999, xhigh = 1.;
    TH1F h_B0_cosAlpha2DwrtBSwithZ_MC = TH1F("B0_cosAlpha2DwrtBSwithZ_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_cosAlpha2DwrtBSwithZ_Fk = TH1F("B0_cosAlpha2DwrtBSwithZ_Fk", "", Nbins, xlow, xhigh);
    TH1F h_B0_cosAlpha2DwrtBS_MC = TH1F("B0_cosAlpha2DwrtBS_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_cosAlpha2DwrtBS_Fk = TH1F("B0_cosAlpha2DwrtBS_Fk", "", Nbins, xlow, xhigh);
    TH1F h_B0_cosAlpha3DwrtPV_MC = TH1F("B0_cosAlpha3DwrtPV_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_cosAlpha3DwrtPV_Fk = TH1F("B0_cosAlpha3DwrtPV_Fk", "", Nbins, xlow, xhigh);
    Nbins = 50, xlow = 0., xhigh = 100;
    TH1F h_B0_LxySign_wrtPV_MC = TH1F("B0_LxySign_wrtPV_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_LxySign_wrtPV_Fk = TH1F("B0_LxySign_wrtPV_Fk", "", Nbins, xlow, xhigh);
    TH1F h_B0_LxySign_wrtBS_MC = TH1F("B0_LxySign_wrtBS_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_LxySign_wrtBS_Fk = TH1F("B0_LxySign_wrtBS_Fk", "", Nbins, xlow, xhigh);
    TH1F h_B0_LxySign_wrtBSwithZ_MC = TH1F("B0_LxySign_wrtBSwithZ_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_LxySign_wrtBSwithZ_Fk = TH1F("B0_LxySign_wrtBSwithZ_Fk", "", Nbins, xlow, xhigh);
    Nbins = 50, xlow = 0., xhigh = 0.015;
    TH1F h_B0_DistGenVtx_PV_MC   = TH1F("B0_DistGenVtx_PV_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_DistGenVtx_PV_Fk   = TH1F("B0_DistGenVtx_PV_Fk", "", Nbins, xlow, xhigh);
    TH1F h_B0_DistGenVtx_BS_MC   = TH1F("B0_DistGenVtx_BS_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_DistGenVtx_BS_Fk   = TH1F("B0_DistGenVtx_BS_Fk", "", Nbins, xlow, xhigh);
    TH1F h_B0_DistGenVtx_BSwithZ_MC = TH1F("B0_DistGenVtx_BSwithZ_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_DistGenVtx_BSwithZ_Fk = TH1F("B0_DistGenVtx_BSwithZ_Fk", "", Nbins, xlow, xhigh);
    Nbins = 20, xlow = 0., xhigh = 1.;
    TH1F h_B0_svProb_MC   = TH1F("B0_svProb_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_svProb_Fk   = TH1F("B0_svProb_Fk", "", Nbins, xlow, xhigh);
    Nbins = 30, xlow = 0., xhigh = 15;
    TH1F h_B0_pT_MC       = TH1F("B0_pT_MC", "", Nbins, xlow, xhigh);
    TH1F h_B0_pT_Fk       = TH1F("B0_pT_Fk", "", Nbins, xlow, xhigh);

    // MCmatching variables
    bool isMCmatched_Mu1, isMCmatched_Mu2, isMCmatched_Pi1, isMCmatched_Pi2;
    bool isMCmatched_JPsi, isMCmatched_PiPi, isMCmatched_Psi2S, isMCmatched_K0s, isMCmatched_B0;
    int N_FiredEvents = 0, N_PassedEvents = 0, n_PassedB0 = 0, N_PassedB0 = 0, N_B0matching = 0;
    // toCount variables
    bool toCountJPsi = true, toCountPiPi= true, toCountK0s = true;
    int  prevMu1_idx, prevMu2_idx, prevPi1_idx, prevPi2_idx;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0 || jentry == Nbreak) break;
        if ((jentry+1) % Nprint == 0) std::cout << "--> " << Form("%3.0f",(float)(jentry+1)/nentries* 100.) << " \%"<< std::endl;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // ----- CHECK IF THE TRIGGER FIRED
        if(!HLT_DoubleMu4_3_LowMass) continue;
        N_FiredEvents++;
        // ----- GENERATOR
        GenPartFillP4();
        realB_idx = GenB0idx();

        // ----- FIND THE MONTE CARLO TRUTH
        //std::cout << " --- EV " << jentry << " with #B0 candidates " << nB0 << std::endl;
        MCtruthMatching();

        n_PassedB0 = 0;
        for (Int_t b = 0; b  < nB0; b++){
            
            if ( !RecoPartFillP4(b)) continue; 
            
            // --> check the MC matching of the chain
            isMCmatched_Mu1   =  (B0_mu1_idx[b] == MCmatch_Mum_Idx) || (B0_mu1_idx[b] == MCmatch_Mup_Idx);
            isMCmatched_Mu2   =  (B0_mu2_idx[b] == MCmatch_Mum_Idx) || (B0_mu2_idx[b] == MCmatch_Mup_Idx);
            isMCmatched_JPsi  =  isMCmatched_Mu1 && isMCmatched_Mu2;
            isMCmatched_Pi1   =  (B0_pi1_idx[b] == MCmatch_Pim_Idx) || (B0_pi1_idx[b] == MCmatch_Pip_Idx);
            isMCmatched_Pi2   =  (B0_pi2_idx[b] == MCmatch_Pim_Idx) || (B0_pi2_idx[b] == MCmatch_Pip_Idx);
            isMCmatched_PiPi   =  isMCmatched_Pi1 && isMCmatched_Pi2;
            isMCmatched_Psi2S =  isMCmatched_PiPi && isMCmatched_JPsi;
            isMCmatched_K0s   =  (fabs(ROOT::Math::VectorUtil::DeltaR(GenP4_K0s, RecoP4_K0s) - MCmatch_K0s_DRmin) < 0.0001) && (MCmatch_K0s_Idx  != -1);
            isMCmatched_B0    =  isMCmatched_Psi2S &&   isMCmatched_K0s;   

            // ----- TRIGGER SELECTION TO B0 CANDIDATE
            if(!TriggerSelection_Muons(b)) continue;
            if(!TriggerSelection_Track(b)) continue;
            n_PassedB0++;

            // *** JPsi --> Mu+Mu-  ***
            toCountJPsi = (B0_mu1_idx[b] != prevMu1_idx) || (B0_mu2_idx[b] != prevMu2_idx);
            if(toCountJPsi){
                prevMu1_idx = B0_mu1_idx[b]; prevMu2_idx = B0_mu2_idx[b];

                if (isMCmatched_Mu1){
                    h_Mu_SoftID_MC.Fill(Muon_softId[B0_mu1_idx[b]]);
                    h_Mu_GlobalMu_MC.Fill(Muon_isGlobal[B0_mu1_idx[b]]);
                    if (Muon_trackQuality[B0_mu1_idx[b]] < 2) h_Mu_TrkQlty_MC.Fill(Muon_trackQuality[B0_mu1_idx[b]]);
                    else h_Mu_TrkQlty_MC.Fill(2);

                    //if (HLT_Dimuon25_Jpsi) h_Mu_dR_HLT_Dimuon25_Jpsi_MC.Fill(B0_MuMu_mu1_dr_Dimuon25_Jpsi[b]);
                    if (HLT_DoubleMu4_3_LowMass) h_Mu_dR_HLT_DoubleMu4_JpsiTrk_MC.Fill(B0_MuMu_mu1_dr_DoubleMu4_JpsiTrk_Displaced[b]);
                }else{
                    h_Mu_SoftID_Fk.Fill(Muon_softId[B0_mu1_idx[b]]);
                    h_Mu_GlobalMu_Fk.Fill(Muon_isGlobal[B0_mu1_idx[b]]);
                    if (Muon_trackQuality[B0_mu1_idx[b]] < 2) h_Mu_TrkQlty_Fk.Fill(Muon_trackQuality[B0_mu1_idx[b]]);
                    else h_Mu_TrkQlty_Fk.Fill(2);

                    //if (HLT_Dimuon25_Jpsi) h_Mu_dR_HLT_Dimuon25_Jpsi_Fk.Fill(B0_MuMu_mu1_dr_Dimuon25_Jpsi[b]);
                    if (HLT_DoubleMu4_3_LowMass) h_Mu_dR_HLT_DoubleMu4_JpsiTrk_Fk.Fill(B0_MuMu_mu1_dr_DoubleMu4_JpsiTrk_Displaced[b]);
                }
                if (isMCmatched_Mu2){
                    h_Mu_SoftID_MC.Fill(Muon_softId[B0_mu2_idx[b]]);
                    h_Mu_GlobalMu_MC.Fill(Muon_isGlobal[B0_mu2_idx[b]]);
                    if (Muon_trackQuality[B0_mu1_idx[b]] < 2) h_Mu_TrkQlty_MC.Fill(Muon_trackQuality[B0_mu2_idx[b]]);
                    else  h_Mu_TrkQlty_MC.Fill(2);

                    //h_Mu_dR_HLT_Dimuon25_Jpsi_MC.Fill(B0_MuMu_mu2_dr_Dimuon25_Jpsi[b]);
                    h_Mu_dR_HLT_DoubleMu4_JpsiTrk_MC.Fill(B0_MuMu_mu2_dr_DoubleMu4_JpsiTrk_Displaced[b]);
                }else{
                    h_Mu_SoftID_Fk.Fill(Muon_softId[B0_mu2_idx[b]]);
                    h_Mu_GlobalMu_Fk.Fill(Muon_isGlobal[B0_mu2_idx[b]]);
                    if (Muon_trackQuality[B0_mu1_idx[b]] < 2) h_Mu_TrkQlty_Fk.Fill(Muon_trackQuality[B0_mu2_idx[b]]);
                    else h_Mu_TrkQlty_Fk.Fill(2);

                    //h_Mu_dR_HLT_Dimuon25_Jpsi_Fk.Fill(B0_MuMu_mu2_dr_Dimuon25_Jpsi[b]);
                    h_Mu_dR_HLT_DoubleMu4_JpsiTrk_Fk.Fill(B0_MuMu_mu2_dr_DoubleMu4_JpsiTrk_Displaced[b]);
                }

                if(isMCmatched_Mu1 && isMCmatched_Mu2) h_MuMu_M_MC.Fill(((RecoP4_Mu1+ RecoP4_Mu2).M()));
                else h_MuMu_M_Fk.Fill(((RecoP4_Mu1+ RecoP4_Mu2).M()));

            }

        
            // *** Rho --> Pi+Pi-  ***
            toCountPiPi = (prevPi1_idx != B0_pi1_idx[b]) || (prevPi2_idx != B0_pi2_idx[b]) ; 
            if(toCountPiPi){
                prevPi1_idx = B0_pi1_idx[b]; prevPi2_idx = B0_pi2_idx[b];

                if(isMCmatched_Pi1){
                    if (HLT_DoubleMu4_3_LowMass) h_Pi_dR_HLT_DoubleMu4_JpsiTrk_MC.Fill(B0_PiPi_p1_dr_DoubleMu4_JpsiTrk_Displaced[b]);
                }else{
                    if (HLT_DoubleMu4_3_LowMass) h_Pi_dR_HLT_DoubleMu4_JpsiTrk_Fk.Fill(B0_PiPi_p1_dr_DoubleMu4_JpsiTrk_Displaced[b]);
                }
                if (isMCmatched_Pi2)
                {
                    if (HLT_DoubleMu4_3_LowMass) h_Pi_dR_HLT_DoubleMu4_JpsiTrk_MC.Fill(B0_PiPi_p2_dr_DoubleMu4_JpsiTrk_Displaced[b]);   
                }else{
                    if (HLT_DoubleMu4_3_LowMass) h_Pi_dR_HLT_DoubleMu4_JpsiTrk_Fk.Fill(B0_PiPi_p2_dr_DoubleMu4_JpsiTrk_Displaced[b]);
                }
                if (isMCmatched_PiPi) h_PiPi_M_MC.Fill( (RecoP4_Pi1+ RecoP4_Pi2).M() ); 
                else h_PiPi_M_Fk.Fill( (RecoP4_Pi1+ RecoP4_Pi2).M() );

            }

        // *** K0short ***
        if (isMCmatched_K0s){
 
        }else{
            

        }
        

        Dist_GenV_PV = sqrt( (GenPart_vx[realB_idx] - B0_PVx[b])*(GenPart_vx[realB_idx] - B0_PVx[b]) + (GenPart_vy[realB_idx] - B0_PVy[b])*(GenPart_vy[realB_idx] - B0_PVy[b]));
        Dist_GenV_BS = sqrt( (GenPart_vx[realB_idx] - B0_BSxRaw[b])*(GenPart_vx[realB_idx] - B0_BSxRaw[b]) + (GenPart_vy[realB_idx] - B0_BSyRaw[b])*(GenPart_vy[realB_idx] - B0_BSyRaw[b]));
        Dist_GenV_BSwithZ = sqrt( (GenPart_vx[realB_idx] - B0_BSxWithZ[b])*(GenPart_vx[realB_idx] - B0_BSxWithZ[b]) + (GenPart_vy[realB_idx] - B0_BSyWithZ[b])*(GenPart_vy[realB_idx] - B0_BSyWithZ[b]));

        // *** B0 --> Psi(3872) K0s ***
        if (isMCmatched_B0){
            
            // ... masses ...
            M_MuMu = B0_MuMu_fitted_mass[b];
            M_PiPi = B0_finalFit_Rho_mass[b];
            M_Psi2S= B0_finalFit_X_mass[b];
            M_K0s = B0_K0s_nmcFitted_mass[b];
            M_B0 = B0_finalFit_mass[b];
            // ... pT ...
            pT_Mu1 = RecoP4_Mu1.Pt(); 
            pT_Mu2 = RecoP4_Mu2.Pt(); 
            pT_Pi1 = RecoP4_Pi1.Pt(); 
            pT_Pi2 = RecoP4_Pi2.Pt(); 
            pT_K0s = RecoP4_K0s.Pt(); 
            
            outTree_->Fill();

            // Rho
            h_Pi1_pT_MC.Fill(RecoP4_Pi1.Pt()/RecoP4_B0.Pt());
            h_Pi1_D0_MC.Fill( B0_PiPi_pi1_d0sig[b] );
            h_Pi1_DRwrtB_MC.Fill(ROOT::Math::VectorUtil::DeltaR(RecoP4_Pi1, RecoP4_B0));
            h_PiPi_svProb_MC.Fill(B0_PiPi_sv_prob[b]);
            h_PiPi_pT_MC.Fill( (RecoP4_Pi1 + RecoP4_Pi2).Pt()/RecoP4_B0.Pt() );

            //K0short
            h_K0s_LxySign_wrtBvtx_MC.Fill(B0_K0_lxySign_wrtBvtx[b]);
            //B0
            h_B0_cosAlpha3DwrtPV_MC.Fill(fabs(B0_cosAlpha3D_PV[b]));
            h_B0_cosAlpha2DwrtBS_MC.Fill(fabs(B0_cosAlpha2D_BS[b]));
            h_B0_cosAlpha2DwrtBSwithZ_MC.Fill(fabs(B0_cosAlpha2D_BSwithZ[b]));
            

            h_B0_LxySign_wrtPV_MC.Fill(B0_lxySign_PV[b]);
            h_B0_LxySign_wrtBS_MC.Fill(B0_lxySign_BS[b]);
            h_B0_LxySign_wrtBSwithZ_MC.Fill(B0_lxySign_BSwithZ[b]);

            h_B0_svProb_MC.Fill(B0_svprob[b]);
            h_B0_pT_MC.Fill(RecoP4_B0.Pt()/B0_finalFit_mass[b]);
           
            h_B0_DistGenVtx_PV_MC.Fill(Dist_GenV_PV);
            h_B0_DistGenVtx_BS_MC.Fill(Dist_GenV_BS);
            h_B0_DistGenVtx_BSwithZ_MC.Fill(Dist_GenV_BSwithZ);
        }else{

            // Rho
            h_Pi1_pT_Fk.Fill(RecoP4_Pi1.Pt()/RecoP4_B0.Pt());
            h_Pi1_D0_Fk.Fill( B0_PiPi_pi1_d0sig[b] );
            h_Pi1_DRwrtB_Fk.Fill(ROOT::Math::VectorUtil::DeltaR(RecoP4_Pi1, RecoP4_B0));
            h_PiPi_svProb_Fk.Fill(B0_PiPi_sv_prob[b]);
            h_PiPi_pT_Fk.Fill( (RecoP4_Pi1 + RecoP4_Pi2).Pt()/RecoP4_B0.Pt() );

            //K0short
            h_K0s_LxySign_wrtBvtx_Fk.Fill(B0_K0_lxySign_wrtBvtx[b]);

            h_B0_cosAlpha3DwrtPV_Fk.Fill(fabs(B0_cosAlpha3D_PV[b]));
            h_B0_cosAlpha2DwrtBS_Fk.Fill(fabs(B0_cosAlpha2D_BS[b]));
            h_B0_cosAlpha2DwrtBSwithZ_Fk.Fill(fabs(B0_cosAlpha2D_BSwithZ[b]));

            h_B0_svProb_Fk.Fill(B0_svprob[b]);
            h_B0_pT_Fk.Fill(RecoP4_B0.Pt()/B0_finalFit_mass[b]);
            

            h_B0_LxySign_wrtPV_Fk.Fill(B0_lxySign_PV[b]);
            h_B0_LxySign_wrtBS_Fk.Fill(B0_lxySign_BS[b]);
            h_B0_LxySign_wrtBSwithZ_Fk.Fill(B0_lxySign_BSwithZ[b]);
            
            h_B0_DistGenVtx_PV_Fk.Fill(Dist_GenV_PV);
            h_B0_DistGenVtx_BS_Fk.Fill(Dist_GenV_BS);
            h_B0_DistGenVtx_BSwithZ_Fk.Fill(Dist_GenV_BSwithZ);
        
        }

        if (isMCmatched_B0) N_B0matching++;

        
        } // loop on B0 candidates
       
        if(n_PassedB0 > 0){
            N_PassedB0 += n_PassedB0;
            N_PassedEvents++;
        } 

    }// loop on events

    std::cout << " Events which fired the HLT "  << N_FiredEvents << std::endl;
    std::cout << " Events which passed the HLT " << N_PassedEvents << std::endl;
    std::cout << " B0 cand. which passed the HLT " << N_PassedB0 << std::endl;
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

    h_MuMu_M_MC.Write();
    h_MuMu_M_Fk.Write();
    h_PiPi_M_MC.Write();
    h_PiPi_M_Fk.Write();

    h_Pi_dR_HLT_DoubleMu4_JpsiTrk_MC.Write();
    h_Pi_dR_HLT_DoubleMu4_JpsiTrk_Fk.Write();
    h_Pi1_pT_MC.Write();
    h_Pi1_pT_Fk.Write();
    h_Pi1_D0_MC.Write();
    h_Pi1_D0_Fk.Write();
    h_Pi1_DRwrtB_MC.Write();
    h_Pi1_DRwrtB_Fk.Write();
    h_PiPi_pT_MC.Write();
    h_PiPi_pT_Fk.Write();
    h_PiPi_svProb_MC.Write();
    h_PiPi_svProb_Fk.Write();

    h_K0s_LxySign_wrtBvtx_MC.Write();
    h_K0s_LxySign_wrtBvtx_Fk.Write();

    h_B0_cosAlpha2DwrtBS_MC.Write();
    h_B0_cosAlpha2DwrtBS_Fk.Write();
    h_B0_cosAlpha2DwrtBSwithZ_MC.Write();
    h_B0_cosAlpha2DwrtBSwithZ_Fk.Write();
    h_B0_cosAlpha3DwrtPV_MC.Write();
    h_B0_cosAlpha3DwrtPV_Fk.Write();

    h_B0_LxySign_wrtPV_MC.Write();
    h_B0_LxySign_wrtPV_Fk.Write();
    h_B0_LxySign_wrtBS_MC.Write();
    h_B0_LxySign_wrtBS_Fk.Write();
    h_B0_LxySign_wrtBSwithZ_MC.Write();
    h_B0_LxySign_wrtBSwithZ_Fk.Write();
    
    h_B0_svProb_MC.Write();
    h_B0_svProb_Fk.Write();
    h_B0_pT_MC.Write();
    h_B0_pT_Fk.Write();
    h_B0_DistGenVtx_PV_MC.Write();
    h_B0_DistGenVtx_BS_MC.Write();
    h_B0_DistGenVtx_BSwithZ_MC.Write();

    outTree_->Write();

    outFile_->Close();
    std::cout << "  ...[OUTPUT] output histograms written on file " << outFilePath_ << std::endl;


}//Loop()


int RecoDecayPsi::RecoPartFillP4(const int Bidx){
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

    RecoP4_B0.SetPt(B0_finalFit_pt[Bidx]); RecoP4_B0.SetEta(B0_finalFit_eta[Bidx]); RecoP4_B0.SetPhi(B0_finalFit_phi[Bidx]); RecoP4_B0.SetM(B0_finalFit_mass[Bidx]);
 
    return TrackQualityCheck;
 
 }//RecoPartFillP4()

void RecoDecayPsi::MCtruthMatching(const bool verbose){

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


void RecoDecayPsi::OutTree_setup(){

    TString TreeName = "RecoDecay_MCmatch";

	
	outTree_ = new TTree( TreeName, TreeName);
	std::cout << " out tree setting up ... " << std::endl;

	outTree_->Branch("run", &Run, "run/F");
	outTree_->Branch("LumiBlock", &LumiBlock, "LumiBlock/F");
	outTree_->Branch("event", &Event, "Event/F");
	
	outTree_->Branch("M_MuMu", &M_MuMu, "M_MuMu/F");
	outTree_->Branch("M_PiPi", &M_PiPi, "M_PiPi/F");
	outTree_->Branch("M_Psi2S", &M_Psi2S, "M_Psi2S/F");
	outTree_->Branch("M_K0s", &M_K0s, "M_K0s/F");
	outTree_->Branch("M_B0", &M_B0, "M_B0/F");

    outTree_->Branch("pT_Mu1", &pT_Mu1, "pT_Mu1/F");
    outTree_->Branch("pT_Mu2", &pT_Mu2, "pT_Mu2/F");
    outTree_->Branch("pT_Pi1", &pT_Pi1, "pT_Pi1/F");
    outTree_->Branch("pT_Pi2", &pT_Pi2, "pT_Pi2/F");
    outTree_->Branch("pT_K0s", &pT_K0s, "pT_K0s/F");

}//OutTree_setup()

float RecoDecayPsi::DeltaPT(ROOT::Math::PtEtaPhiMVector& genV, ROOT::Math::PtEtaPhiMVector& recV){
	return fabs(genV.Pt() - recV.Pt()) / genV.Pt();
}//DeltaPT()

int RecoDecayPsi::TriggerSelection_Muons(const int Bidx){
   // TRIGGER SETTINGS 
    const float Min_Mu_pT = 4.,Max_Mu_eta = 2.5, Max_Mu_dr = 2.;
    const float Min_MuMu_pT = 6.9, Low_MuMu_M = 3.0,  High_MuMu_M = 3.2, Max_MuMu_DCA = 0.5;
    const float Min_MuMu_LxyS = 3, Min_MuMu_cosAlpha = 0.9, Min_MuMu_SVp = 0.1;

    int mu1_idx, mu2_idx;
    bool isFiredMu1, isFiredMu2;
    bool isOK_mu1_step0 = false, isOK_mu2_step0 = false, MassCut = false, isOK_mumu_step1 = false, isOK_mumu_step2 = false; 

    int RETURN_VALUE = 0;

    mu1_idx = B0_mu1_idx[Bidx];
    mu2_idx = B0_mu2_idx[Bidx];

    // Fired Mu + muon tracks QUALITY CHECK
    isFiredMu1 = (bool)B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];
    isFiredMu2 = (bool)B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced[Bidx]; 
    if ( (isFiredMu1 && isFiredMu2) && ( Muon_softId[mu1_idx] && Muon_softId[mu2_idx] )){ 
            // STEP 0
            isOK_mu1_step0 = true;
            if((RecoP4_Mu1.Pt() < Min_Mu_pT) || (fabs(RecoP4_Mu1.Eta()) > Max_Mu_eta)  || ( B0_MuMu_mu1_dr[Bidx]) > Max_Mu_dr ) isOK_mu1_step0 = false;
            isOK_mu2_step0 = true;
            if((RecoP4_Mu2.Pt() < Min_Mu_pT) || (fabs(RecoP4_Mu2.Eta()) > Max_Mu_eta)  || ( B0_MuMu_mu2_dr[Bidx]) > Max_Mu_dr ) isOK_mu2_step0 = false;

            if (isOK_mu1_step0 && isOK_mu2_step0){ 

                // STEP 1 
                isOK_mumu_step1 = true;
                MassCut = ( (RecoP4_Mu1 + RecoP4_Mu2).M() > Low_MuMu_M ) && ( (RecoP4_Mu1 + RecoP4_Mu2).M() < High_MuMu_M  );
                if ( !MassCut || ((RecoP4_Mu1 + RecoP4_Mu2).Pt() < Min_MuMu_pT ) || ( B0_MuMu_DCA[Bidx] > Max_MuMu_DCA )  )	isOK_mumu_step1 = false;
                // STEP 2	
                isOK_mumu_step2 = true;
                if((B0_MuMu_LxySign[Bidx] < Min_MuMu_LxyS) || (B0_MuMu_cosAlpha[Bidx] < Min_MuMu_cosAlpha ) || (B0_MuMu_sv_prob[Bidx] < Min_MuMu_SVp )) isOK_mumu_step2 = false;
            }
    }

    if (isOK_mu1_step0 && isOK_mu2_step0 && isOK_mumu_step1 && isOK_mumu_step2) RETURN_VALUE = 1; 

    return RETURN_VALUE;

}//TriggerSelection_Muons()

int RecoDecayPsi::TriggerSelection_Track(const int Bidx){
   //TRIGGER SETTINGS
   const float Min_Trk_pT = 1.2, Max_Trk_eta = 2.5, Min_Trk_D0S = 2.;
   bool isOK_trk_step0 = false, isOK_trk_step1 = false;

   int RETURN_VALUE = 0;

	bool isFired_RhoPi1 = (bool)B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];
	bool isMatchedToMuon_Rho_Pi1	= ProbeTracks_isMatchedToMuon[B0_pi1_idx[Bidx]];
	bool isFired_RhoPi2 = (bool)B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];
	bool isMatchedToMuon_Rho_Pi2 = ProbeTracks_isMatchedToMuon[B0_pi2_idx[Bidx]];

	bool isFired_K0sPi1 = (bool)B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];
	bool isFired_K0sPi2 = (bool)B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];

	//LEVEL 0 
	isOK_trk_step0 =  (isFired_RhoPi1 || isFired_RhoPi2 || isFired_K0sPi1 || isFired_K0sPi2); // is fired at least 1

   // LEVEL 1
   isOK_trk_step1 = false;
   if (isFired_RhoPi1 && !isMatchedToMuon_Rho_Pi1){ 
		if( (RecoP4_Pi1.Pt() > Min_Trk_pT) && (fabs(RecoP4_Pi1.Eta()) < Max_Trk_eta) && (B0_PiPi_pi1_d0sig[Bidx] > Min_Trk_D0S) ){
			isOK_trk_step1 = true;
			RETURN_VALUE = 1;
		}
   } 
   if (isFired_RhoPi2 && !isMatchedToMuon_Rho_Pi2){ 
      if( (RecoP4_Pi2.Pt() > Min_Trk_pT) && (fabs(RecoP4_Pi2.Eta()) < Max_Trk_eta) && (B0_PiPi_pi2_d0sig[Bidx] > Min_Trk_D0S) ){
			isOK_trk_step1 = true;
			RETURN_VALUE = 2;
		}
   }
   if (isFired_K0sPi1){
      if( (B0_K0s_matchTrack1_pt[Bidx] > Min_Trk_pT) && (fabs(B0_K0s_matchTrack1_eta[Bidx]) < Max_Trk_eta) && (B0_K0s_matchTrack1_D0sign[Bidx] > Min_Trk_D0S) ){
			isOK_trk_step1 = true;
			RETURN_VALUE = 3;
		}
   }
   if (isFired_K0sPi2){
      if( (B0_K0s_matchTrack2_pt[Bidx] > Min_Trk_pT) && (fabs(B0_K0s_matchTrack2_eta[Bidx]) < Max_Trk_eta) && (B0_K0s_matchTrack2_D0sign[Bidx] > Min_Trk_D0S) ){
			isOK_trk_step1 = true;
			RETURN_VALUE = 3;
		}
   }

	//if (isOK_trk_step0 && isOK_trk_step1) RETURN_VALUE = 1;	

	return RETURN_VALUE;

}//TriggerSelection_Tracks
