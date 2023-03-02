#define B0toX3872K0s_base_cxx
#include "../include/B0toX3872K0s_base.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void B0toX3872K0s_base::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

int B0toX3872K0s_base::RecoPartFillP4(const int Bidx){

   int TrackQualityCheck = 1;

   RecoP4_Mu1.SetM(mMuon); RecoP4_Mu2.SetM(mMuon);
   RecoP4_Pi1.SetM(mPion); RecoP4_Pi2.SetM(mPion);
   RecoP4_K0s.SetM(mK0s);

   //... muons P4
   if(!Muon_softId[B0_mu1_idx[Bidx]] || !Muon_softId[B0_mu2_idx[Bidx]]) TrackQualityCheck = 0;
   RecoP4_Mu1.SetPt(B0_finalFit_mu1_pt[Bidx]); RecoP4_Mu1.SetEta(B0_finalFit_mu1_eta[Bidx]); RecoP4_Mu1.SetPhi(B0_finalFit_mu1_phi[Bidx]);
   RecoP4_Mu2.SetPt(B0_finalFit_mu2_pt[Bidx]); RecoP4_Mu2.SetEta(B0_finalFit_mu2_eta[Bidx]); RecoP4_Mu2.SetPhi(B0_finalFit_mu2_phi[Bidx]);
   
   RecoP4_JPsi = RecoP4_Mu1 + RecoP4_Mu2;
   
   //... tracks P4
   if(ProbeTracks_isMatchedToMuon[B0_pi1_idx[Bidx]] || ProbeTracks_isMatchedToMuon[B0_pi2_idx[Bidx]])TrackQualityCheck = 0;
   RecoP4_Pi1.SetPt(B0_finalFit_pi1_pt[Bidx]); RecoP4_Pi1.SetEta(B0_finalFit_pi1_eta[Bidx]); RecoP4_Pi1.SetPhi(B0_finalFit_pi1_phi[Bidx]);
   RecoP4_Pi2.SetPt(B0_finalFit_pi2_pt[Bidx]); RecoP4_Pi2.SetEta(B0_finalFit_pi2_eta[Bidx]); RecoP4_Pi2.SetPhi(B0_finalFit_pi2_phi[Bidx]);
   
   RecoP4_PiPi = RecoP4_Pi1 + RecoP4_Pi2;
   RecoP4_PiPi.SetM(B0_finalFit_Rho_mass[Bidx]);
   RecoP4_X3872 = RecoP4_JPsi + RecoP4_PiPi;
   RecoP4_X3872.SetM(B0_finalFit_Rho_mass[Bidx]);
   
   RecoP4_K0s.SetPt(B0_K0s_mcFitted_pt[Bidx]); RecoP4_K0s.SetEta(B0_K0s_mcFitted_eta[Bidx]); RecoP4_K0s.SetPhi(B0_K0s_mcFitted_phi[Bidx]);

   RecoP4_B0.SetPt(B0_finalFit_pt[Bidx]); RecoP4_B0.SetEta(B0_finalFit_eta[Bidx]); RecoP4_B0.SetPhi(B0_finalFit_phi[Bidx]); RecoP4_B0.SetM(B0_finalFit_mass[Bidx]);

   return TrackQualityCheck;

}