#include "../include/HLTapply.h"


HLTapply::HLTapply(TTree *tree, const TString channel, const TString tags) : B0toX3872K0s_base(tree, channel) {
    tags_ = tags;
    TString blind_tag = "Blind";
    if(!isBlind_) blind_tag = "Open";
    outFileTreePath_ = "outRoot/Data" + blind_tag + "_HLTemulation_" + channel_ + "_" + tags_ + ".root";
    
    //std::cout << " .... analyzing " << tree->GetEntriesFast() << " events "<< std::endl;

}

void HLTapply::Loop(){

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    Long64_t Nbreak =  nentries + 10, Nprint = 1000;//(int)((double)nentries/20.); 
    long int TotalEvents = 0;


    // ----- VARIABLES ----- //
    // toCount variables
    int  N_FiredEvents = 0, N_PassedEvents = 0, n_PassedB0 = 0, N_PassedB0 = 0, N_B0matching = 0;
    bool toCountJPsi = true, toCountPiPi= true, toCountK0s = true;
    int  prevMu1_idx, prevMu2_idx, prevPi1_idx, prevPi2_idx;

    ROOT::Math::PtEtaPhiMVector prevRecoP4_K0s(0,0,0,0);

    // ----- OUTPUT TREE SETUP ----- //
	OutTree_setup();

    for (Long64_t jentry=0; jentry<nentries;jentry++) {

        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0 || jentry == Nbreak) break;
        TotalEvents++;

        if ((jentry+1) % Nprint == 0) std::cout << " --> " <<  (jentry+1) << std::endl;//std::cout << "--> " << Form("%3.0f",(float)(jentry+1)/nentries* 100.) << " \%"<< std::endl;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        
        // ----- CHECK IF THE TRIGGER FIRED
        if(!HLT_DoubleMu4_JpsiTrk_Displaced) continue;
        N_FiredEvents++;

        for (Int_t b = 0; b  < nB0; b++){
            
            // blind the B0 mass region if specified
            M_B0 = B0_finalFit_mass[b];
            M_X3872 = B0_finalFit_X_mass[b];
            if (isBlind_ && (M_B0 > Blind_MB0_low && M_B0 < Blind_MB0_high && M_X3872 > 3.75) ) continue;

            // save the reconstructed momenta for the particle and check te quality tracks
            if ( !RecoPartFillP4(b)) continue;


            // ----- TRIGGER SELECTION TO THE B0 CANDIDATE
            //std::cout << " - muon selection " << TriggerSelection_Muons(b) << std::endl;
            //std::cout << " - track selection " << TriggerSelection_Track(b) << std::endl;
            if(!TriggerSelection_Muons(b)) continue;
            if(!TriggerSelection_Track(b)) continue;
            n_PassedB0++;

            toCountJPsi = (B0_mu1_idx[b] != prevMu1_idx) || (B0_mu2_idx[b] != prevMu2_idx);
            if (toCountJPsi){
                prevMu1_idx = B0_mu1_idx[b]; prevMu2_idx = B0_mu2_idx[b];
            }

            toCountPiPi = (B0_pi1_idx[b] != prevPi1_idx) || (B0_pi2_idx[b] != prevPi2_idx);
            if (toCountPiPi){
                prevPi1_idx = B0_pi1_idx[b]; prevPi2_idx = B0_pi2_idx[b];
            }

            toCountK0s = ROOT::Math::VectorUtil::DeltaR(prevRecoP4_K0s, RecoP4_K0s)  < 0.0001;
            if (toCountK0s) prevRecoP4_K0s = RecoP4_K0s;


        // ********** output treee **********
            Run = run;
            LumiBlock = luminosityBlock;
            Event = event;

            M_mumu = B0_MuMu_fitted_mass[b];
            M_PiPi = B0_finalFit_Rho_mass[b];
            M_X3872 = B0_finalFit_X_mass[b];
            M_K0sPi1 = RecoP4_K0sPi1.M();
            M_K0sPi2 = RecoP4_K0sPi2.M();
            M_K0s = B0_K0s_nmcFitted_mass[b];


            			
            LxySignBSz_B0       = B0_lxySign_BSwithZ[b];
            SVprob_B0           = B0_svprob[b];
            CosAlpha3DBSz_B0    = B0_cosAlpha2D_BSwithZ[b];
            pTM_B0              = RecoP4_B0.Pt()/M_B0;
            LxySignSV_K0s       = B0_K0_lxySign_wrtBvtx[b];
            SVprob_PiPi         = B0_PiPi_sv_prob[b];
            pT_PiPi             = RecoP4_PiPi.Pt()/RecoP4_B0.Pt();
            pT_Pi1              = RecoP4_Pi1.Pt()/RecoP4_B0.Pt();
            DR_B0Pi1            = ROOT::Math::VectorUtil::DeltaR(RecoP4_Pi1, RecoP4_B0);
            D0_Pi1              = B0_PiPi_pi1_d0sig[b];
            
            outTree_->Fill();
            
            

        } // loop on B0 candiadates
        

        if(n_PassedB0 > 0){
            N_PassedB0 += n_PassedB0;
            N_PassedEvents++;
        } 
        n_PassedB0 = 0;

    }// loop on events

    std::cout << " Total processed events "  << TotalEvents << std::endl;
    std::cout << " Events which fired the HLT "  << N_FiredEvents << std::endl;
    std::cout << " Events which passed the HLT " << N_PassedEvents << std::endl;
    std::cout << " B0 cand. which passed the HLT " << N_PassedB0 << std::endl;


    outFileTree_ = new TFile( outFileTreePath_, "RECREATE");	
	if (!outFileTree_->IsOpen()) std::cout << "	ERROR: cannot open Histo out-file " << outFileTreePath_ << std::endl;
	else std::cout << " ... [OUTPUT]  " << outFileTreePath_ << std::endl;

	outFileTree_->cd();
	outTree_->Write();
	outFileTree_->Close();

}


void HLTapply::OutTree_setup(){

    TString TreeName = "HLTemulation";

	
	outTree_ = new TTree( TreeName, TreeName);
	std::cout << " out tree setting up ... " << std::endl;

	outTree_->Branch("run", &Run, "run/F");
	outTree_->Branch("LumiBlock", &LumiBlock, "LumiBlock/F");
	outTree_->Branch("event", &Event, "Event/F");
	
	outTree_->Branch("M_B0", &M_B0, "M_B0/F");
	outTree_->Branch("M_PiPi", &M_PiPi, "M_PiPi/F");
	outTree_->Branch("M_X3872", &M_X3872, "M_X3872/F");
	outTree_->Branch("M_K0s", &M_K0s, "M_K0s/F");
	outTree_->Branch("M_K0sPi1", &M_K0sPi1, "M_K0sPi1/F");
	outTree_->Branch("M_K0sPi2", &M_K0sPi2, "M_K0sPi2/F");
	outTree_->Branch("M_mumu", &M_mumu, "M_mumu/F");

	outTree_->Branch("pTM_B0", &pTM_B0, "pTM_B0/F");
	outTree_->Branch("LxySignBSz_B0", &LxySignBSz_B0, "LxySignBSz_B0/F");
	outTree_->Branch("SVprob_B0", &SVprob_B0, "SVprob_B0/F");
	outTree_->Branch("CosAlpha3DBSz_B0", &CosAlpha3DBSz_B0, "CosAlpha3DBSz_B0/F");
    outTree_->Branch("LxySignSV_K0s", &LxySignSV_K0s, "LxySignSV_K0s/F");
    outTree_->Branch("SVprob_PiPi", &SVprob_PiPi, "SVprob_PiPi/F");
	outTree_->Branch("pT_PiPi", &pT_PiPi, "pT_PiPi/F");
	outTree_->Branch("pT_Pi1", &pT_Pi1, "pT_Pi1/F");
	outTree_->Branch("DR_B0Pi1", &DR_B0Pi1, "DR_B0Pi1/F");
	outTree_->Branch("D0_Pi1", &D0_Pi1, "D0_Pi1/F");


}//OutTree_setup()

int HLTapply::TriggerSelection_Muons(const int Bidx){
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

int HLTapply::TriggerSelection_Track(const int Bidx){
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
