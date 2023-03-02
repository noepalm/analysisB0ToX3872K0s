//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 24 12:23:37 2023 by ROOT version 6.14/09
// from TTree Events/Events
// found on file: data/xNANO_data_2023Feb06_135.root
//////////////////////////////////////////////////////////

#ifndef B0toX3872K0s_base_h
#define B0toX3872K0s_base_h

#include <string>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include <TMath.h>
#include "Math/Vector4D.h"
#include <Math/GenVector/VectorUtil.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include "TLorentzVector.h"

// Header file for the classes stored in the TTree if any.

class B0toX3872K0s_base {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          nB0;
   Float_t         B0_BSxRaw[100];   //[nB0]
   Float_t         B0_BSxWithZ[100];   //[nB0]
   Float_t         B0_BSyRaw[100];   //[nB0]
   Float_t         B0_BSyWithZ[100];   //[nB0]
   Float_t         B0_K0_cosAlpha3D[100];   //[nB0]
   Float_t         B0_K0_lxySign_wrtBvtx[100];   //[nB0]
   Float_t         B0_K0s_matchTrack1_D0sign[100];   //[nB0]
   Float_t         B0_K0s_matchTrack1_dR[100];   //[nB0]
   Float_t         B0_K0s_matchTrack1_eta[100];   //[nB0]
   Float_t         B0_K0s_matchTrack1_maxD0Pv[100];   //[nB0]
   Float_t         B0_K0s_matchTrack1_minD0Pv[100];   //[nB0]
   Float_t         B0_K0s_matchTrack1_phi[100];   //[nB0]
   Float_t         B0_K0s_matchTrack1_pt[100];   //[nB0]
   Float_t         B0_K0s_matchTrack2_D0sign[100];   //[nB0]
   Float_t         B0_K0s_matchTrack2_dR[100];   //[nB0]
   Float_t         B0_K0s_matchTrack2_eta[100];   //[nB0]
   Float_t         B0_K0s_matchTrack2_maxD0Pv[100];   //[nB0]
   Float_t         B0_K0s_matchTrack2_minD0Pv[100];   //[nB0]
   Float_t         B0_K0s_matchTrack2_phi[100];   //[nB0]
   Float_t         B0_K0s_matchTrack2_pt[100];   //[nB0]
   Float_t         B0_K0s_mcFitted_eta[100];   //[nB0]
   Float_t         B0_K0s_mcFitted_mass[100];   //[nB0]
   Float_t         B0_K0s_mcFitted_phi[100];   //[nB0]
   Float_t         B0_K0s_mcFitted_pt[100];   //[nB0]
   Float_t         B0_K0s_mcFitted_svprob[100];   //[nB0]
   Float_t         B0_K0s_mcFitted_vtxX[100];   //[nB0]
   Float_t         B0_K0s_mcFitted_vtxXE[100];   //[nB0]
   Float_t         B0_K0s_mcFitted_vtxY[100];   //[nB0]
   Float_t         B0_K0s_mcFitted_vtxYE[100];   //[nB0]
   Float_t         B0_K0s_mcFitted_vtxZ[100];   //[nB0]
   Float_t         B0_K0s_mcFitted_vtxZE[100];   //[nB0]
   Float_t         B0_K0s_nmcFitted_mass[100];   //[nB0]
   Float_t         B0_K0s_nmcFitted_pi1eta[100];   //[nB0]
   Float_t         B0_K0s_nmcFitted_pi1phi[100];   //[nB0]
   Float_t         B0_K0s_nmcFitted_pi1pt[100];   //[nB0]
   Float_t         B0_K0s_nmcFitted_pi2eta[100];   //[nB0]
   Float_t         B0_K0s_nmcFitted_pi2phi[100];   //[nB0]
   Float_t         B0_K0s_nmcFitted_pi2pt[100];   //[nB0]
   Float_t         B0_K0s_prefit_mass[100];   //[nB0]
   Float_t         B0_MuMu_DCA[100];   //[nB0]
   Float_t         B0_MuMu_LxySign[100];   //[nB0]
   Float_t         B0_MuMu_cosAlpha[100];   //[nB0]
   Float_t         B0_MuMu_fitted_eta[100];   //[nB0]
   Float_t         B0_MuMu_fitted_mass[100];   //[nB0]
   Float_t         B0_MuMu_fitted_phi[100];   //[nB0]
   Float_t         B0_MuMu_fitted_pt[100];   //[nB0]
   Float_t         B0_MuMu_fitted_vtxX[100];   //[nB0]
   Float_t         B0_MuMu_fitted_vtxXE[100];   //[nB0]
   Float_t         B0_MuMu_fitted_vtxY[100];   //[nB0]
   Float_t         B0_MuMu_fitted_vtxYE[100];   //[nB0]
   Float_t         B0_MuMu_fitted_vtxZ[100];   //[nB0]
   Float_t         B0_MuMu_fitted_vtxZE[100];   //[nB0]
   Float_t         B0_MuMu_mu1_dr[100];   //[nB0]
   Float_t         B0_MuMu_mu1_dr_Dimuon20_Jpsi[100];   //[nB0]
   Float_t         B0_MuMu_mu1_dr_DoubleMu4_JpsiTrk_Displaced[100];   //[nB0]
   Float_t         B0_MuMu_mu1_dxysign[100];   //[nB0]
   Float_t         B0_MuMu_mu1_dzsign[100];   //[nB0]
   Float_t         B0_MuMu_mu2_dr[100];   //[nB0]
   Float_t         B0_MuMu_mu2_dr_Dimuon20_Jpsi[100];   //[nB0]
   Float_t         B0_MuMu_mu2_dr_DoubleMu4_JpsiTrk_Displaced[100];   //[nB0]
   Float_t         B0_MuMu_mu2_dxysign[100];   //[nB0]
   Float_t         B0_MuMu_mu2_dzsign[100];   //[nB0]
   Float_t         B0_MuMu_prefit_mu1_eta[100];   //[nB0]
   Float_t         B0_MuMu_prefit_mu1_phi[100];   //[nB0]
   Float_t         B0_MuMu_prefit_mu1_pt[100];   //[nB0]
   Float_t         B0_MuMu_prefit_mu2_eta[100];   //[nB0]
   Float_t         B0_MuMu_prefit_mu2_phi[100];   //[nB0]
   Float_t         B0_MuMu_prefit_mu2_pt[100];   //[nB0]
   Float_t         B0_MuMu_sv_prob[100];   //[nB0]
   Float_t         B0_PVEx[100];   //[nB0]
   Float_t         B0_PVEy[100];   //[nB0]
   Float_t         B0_PVEz[100];   //[nB0]
   Float_t         B0_PVx[100];   //[nB0]
   Float_t         B0_PVy[100];   //[nB0]
   Float_t         B0_PVz[100];   //[nB0]
   Float_t         B0_PiPi_p1_dr_DoubleMu4_JpsiTrk_Displaced[100];   //[nB0]
   Float_t         B0_PiPi_p2_dr_DoubleMu4_JpsiTrk_Displaced[100];   //[nB0]
   Float_t         B0_PiPi_pi1_d0sig[100];   //[nB0]
   Float_t         B0_PiPi_pi1_dxysign[100];   //[nB0]
   Float_t         B0_PiPi_pi1_dzsign[100];   //[nB0]
   Float_t         B0_PiPi_pi1_maxd0PV[100];   //[nB0]
   Float_t         B0_PiPi_pi1_mind0PV[100];   //[nB0]
   Float_t         B0_PiPi_pi2_d0sig[100];   //[nB0]
   Float_t         B0_PiPi_pi2_dxysign[100];   //[nB0]
   Float_t         B0_PiPi_pi2_dzsign[100];   //[nB0]
   Float_t         B0_PiPi_pi2_maxd0PV[100];   //[nB0]
   Float_t         B0_PiPi_pi2_mind0PV[100];   //[nB0]
   Float_t         B0_PiPi_prefit_pi1_eta[100];   //[nB0]
   Float_t         B0_PiPi_prefit_pi1_phi[100];   //[nB0]
   Float_t         B0_PiPi_prefit_pi1_pt[100];   //[nB0]
   Float_t         B0_PiPi_prefit_pi1_vx[100];   //[nB0]
   Float_t         B0_PiPi_prefit_pi1_vy[100];   //[nB0]
   Float_t         B0_PiPi_prefit_pi1_vz[100];   //[nB0]
   Float_t         B0_PiPi_prefit_pi2_eta[100];   //[nB0]
   Float_t         B0_PiPi_prefit_pi2_phi[100];   //[nB0]
   Float_t         B0_PiPi_prefit_pi2_pt[100];   //[nB0]
   Float_t         B0_PiPi_prefit_pi2_vx[100];   //[nB0]
   Float_t         B0_PiPi_prefit_pi2_vy[100];   //[nB0]
   Float_t         B0_PiPi_prefit_pi2_vz[100];   //[nB0]
   Float_t         B0_PiPi_sv_prob[100];   //[nB0]
   Float_t         B0_cosAlpha2D_BS[100];   //[nB0]
   Float_t         B0_cosAlpha2D_BSwithZ[100];   //[nB0]
   Float_t         B0_cosAlpha3D_PV[100];   //[nB0]
   Float_t         B0_decayVtxX[100];   //[nB0]
   Float_t         B0_decayVtxXE[100];   //[nB0]
   Float_t         B0_decayVtxY[100];   //[nB0]
   Float_t         B0_decayVtxYE[100];   //[nB0]
   Float_t         B0_decayVtxZ[100];   //[nB0]
   Float_t         B0_decayVtxZE[100];   //[nB0]
   Float_t         B0_finalFit_JPsi_mass[100];   //[nB0]
   Float_t         B0_finalFit_Rho_mass[100];   //[nB0]
   Float_t         B0_finalFit_X_mass[100];   //[nB0]
   Float_t         B0_finalFit_eta[100];   //[nB0]
   Float_t         B0_finalFit_k0s_eta[100];   //[nB0]
   Float_t         B0_finalFit_k0s_phi[100];   //[nB0]
   Float_t         B0_finalFit_k0s_pt[100];   //[nB0]
   Float_t         B0_finalFit_mass[100];   //[nB0]
   Float_t         B0_finalFit_mu1_eta[100];   //[nB0]
   Float_t         B0_finalFit_mu1_phi[100];   //[nB0]
   Float_t         B0_finalFit_mu1_pt[100];   //[nB0]
   Float_t         B0_finalFit_mu2_eta[100];   //[nB0]
   Float_t         B0_finalFit_mu2_phi[100];   //[nB0]
   Float_t         B0_finalFit_mu2_pt[100];   //[nB0]
   Float_t         B0_finalFit_phi[100];   //[nB0]
   Float_t         B0_finalFit_pi1_eta[100];   //[nB0]
   Float_t         B0_finalFit_pi1_phi[100];   //[nB0]
   Float_t         B0_finalFit_pi1_pt[100];   //[nB0]
   Float_t         B0_finalFit_pi2_eta[100];   //[nB0]
   Float_t         B0_finalFit_pi2_phi[100];   //[nB0]
   Float_t         B0_finalFit_pi2_pt[100];   //[nB0]
   Float_t         B0_finalFit_pt[100];   //[nB0]
   Float_t         B0_fitted_mass_womc[100];   //[nB0]
   Float_t         B0_lxySign_BS[100];   //[nB0]
   Float_t         B0_lxySign_BSwithZ[100];   //[nB0]
   Float_t         B0_lxySign_PV[100];   //[nB0]
   Float_t         B0_svchi2[100];   //[nB0]
   Float_t         B0_svprob[100];   //[nB0]
   Int_t           B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced[100];   //[nB0]
   Int_t           B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced[100];   //[nB0]
   Int_t           B0_MuMu_mu1_fired_Dimuon20_Jpsi[100];   //[nB0]
   Int_t           B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced[100];   //[nB0]
   Int_t           B0_MuMu_mu1_trackQuality[100];   //[nB0]
   Int_t           B0_MuMu_mu2_fired_Dimuon20_Jpsi[100];   //[nB0]
   Int_t           B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced[100];   //[nB0]
   Int_t           B0_MuMu_mu2_trackQuality[100];   //[nB0]
   Int_t           B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced[100];   //[nB0]
   Int_t           B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced[100];   //[nB0]
   Int_t           B0_dimuon_idx[100];   //[nB0]
   Int_t           B0_dipion_idx[100];   //[nB0]
   Int_t           B0_k0short_idx[100];   //[nB0]
   Int_t           B0_mu1_idx[100];   //[nB0]
   Int_t           B0_mu2_idx[100];   //[nB0]
   Int_t           B0_pi1_idx[100];   //[nB0]
   Int_t           B0_pi2_idx[100];   //[nB0]
   Int_t           B0_pv3D_idx[100];   //[nB0]
   UInt_t          nJPsiToMuMu;
   UInt_t          nK0s;
   UInt_t          npipi;
   UInt_t          nMuon;
   Int_t           Muon_trackQuality[100];   //[nMuon]
   Bool_t          Muon_charge[100];   //[nMuon]
   Bool_t          Muon_isGlobal[100];   //[nMuon]
   Bool_t          Muon_looseId[100];   //[nMuon]
   Bool_t          Muon_softId[100];   //[nMuon]
   Float_t         fixedGridRhoFastjetAll;
   Float_t         fixedGridRhoFastjetCentral;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   UInt_t          nProbeTracks;
   Int_t           ProbeTracks_nValidHits[1500];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToLooseMuon[1500];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToMuon[1500];   //[nProbeTracks]
   Bool_t          ProbeTracks_isMatchedToSoftMuon[1500];   //[nProbeTracks]
   Int_t           HLT_Dimuon25_Jpsi;
   Int_t           HLT_DoubleMu4_JpsiTrk_Displaced;
   UInt_t          nTrigObj;
   Float_t         TrigObj_pt[1];   //[nTrigObj]
   Float_t         TrigObj_eta[1];   //[nTrigObj]
   Float_t         TrigObj_phi[1];   //[nTrigObj]
   Int_t           TrigObj_id[1];   //[nTrigObj]
   UInt_t          nOtherPV;
   Float_t         OtherPV_z[10];   //[nOtherPV]
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   UInt_t          nSV;
   Float_t         SV_dlen[10];   //[nSV]
   Float_t         SV_dlenSig[10];   //[nSV]
   Float_t         SV_dxy[10];   //[nSV]
   Float_t         SV_dxySig[10];   //[nSV]
   Float_t         SV_pAngle[10];   //[nSV]
   Int_t           SV_charge[10];   //[nSV]
   Float_t         SV_chi2[10];   //[nSV]
   Float_t         SV_eta[10];   //[nSV]
   Float_t         SV_mass[10];   //[nSV]
   Float_t         SV_ndof[10];   //[nSV]
   Float_t         SV_phi[10];   //[nSV]
   Float_t         SV_pt[10];   //[nSV]
   Float_t         SV_x[10];   //[nSV]
   Float_t         SV_y[10];   //[nSV]
   Float_t         SV_z[10];   //[nSV]
   UChar_t         SV_ntracks[10];   //[nSV]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nB0;   //!
   TBranch        *b_B0_BSxRaw;   //!
   TBranch        *b_B0_BSxWithZ;   //!
   TBranch        *b_B0_BSyRaw;   //!
   TBranch        *b_B0_BSyWithZ;   //!
   TBranch        *b_B0_K0_cosAlpha3D;   //!
   TBranch        *b_B0_K0_lxySign_wrtBvtx;   //!
   TBranch        *b_B0_K0s_matchTrack1_D0sign;   //!
   TBranch        *b_B0_K0s_matchTrack1_dR;   //!
   TBranch        *b_B0_K0s_matchTrack1_eta;   //!
   TBranch        *b_B0_K0s_matchTrack1_maxD0Pv;   //!
   TBranch        *b_B0_K0s_matchTrack1_minD0Pv;   //!
   TBranch        *b_B0_K0s_matchTrack1_phi;   //!
   TBranch        *b_B0_K0s_matchTrack1_pt;   //!
   TBranch        *b_B0_K0s_matchTrack2_D0sign;   //!
   TBranch        *b_B0_K0s_matchTrack2_dR;   //!
   TBranch        *b_B0_K0s_matchTrack2_eta;   //!
   TBranch        *b_B0_K0s_matchTrack2_maxD0Pv;   //!
   TBranch        *b_B0_K0s_matchTrack2_minD0Pv;   //!
   TBranch        *b_B0_K0s_matchTrack2_phi;   //!
   TBranch        *b_B0_K0s_matchTrack2_pt;   //!
   TBranch        *b_B0_K0s_mcFitted_eta;   //!
   TBranch        *b_B0_K0s_mcFitted_mass;   //!
   TBranch        *b_B0_K0s_mcFitted_phi;   //!
   TBranch        *b_B0_K0s_mcFitted_pt;   //!
   TBranch        *b_B0_K0s_mcFitted_svprob;   //!
   TBranch        *b_B0_K0s_mcFitted_vtxX;   //!
   TBranch        *b_B0_K0s_mcFitted_vtxXE;   //!
   TBranch        *b_B0_K0s_mcFitted_vtxY;   //!
   TBranch        *b_B0_K0s_mcFitted_vtxYE;   //!
   TBranch        *b_B0_K0s_mcFitted_vtxZ;   //!
   TBranch        *b_B0_K0s_mcFitted_vtxZE;   //!
   TBranch        *b_B0_K0s_nmcFitted_mass;   //!
   TBranch        *b_B0_K0s_nmcFitted_pi1eta;   //!
   TBranch        *b_B0_K0s_nmcFitted_pi1phi;   //!
   TBranch        *b_B0_K0s_nmcFitted_pi1pt;   //!
   TBranch        *b_B0_K0s_nmcFitted_pi2eta;   //!
   TBranch        *b_B0_K0s_nmcFitted_pi2phi;   //!
   TBranch        *b_B0_K0s_nmcFitted_pi2pt;   //!
   TBranch        *b_B0_K0s_prefit_mass;   //!
   TBranch        *b_B0_MuMu_DCA;   //!
   TBranch        *b_B0_MuMu_LxySign;   //!
   TBranch        *b_B0_MuMu_cosAlpha;   //!
   TBranch        *b_B0_MuMu_fitted_eta;   //!
   TBranch        *b_B0_MuMu_fitted_mass;   //!
   TBranch        *b_B0_MuMu_fitted_phi;   //!
   TBranch        *b_B0_MuMu_fitted_pt;   //!
   TBranch        *b_B0_MuMu_fitted_vtxX;   //!
   TBranch        *b_B0_MuMu_fitted_vtxXE;   //!
   TBranch        *b_B0_MuMu_fitted_vtxY;   //!
   TBranch        *b_B0_MuMu_fitted_vtxYE;   //!
   TBranch        *b_B0_MuMu_fitted_vtxZ;   //!
   TBranch        *b_B0_MuMu_fitted_vtxZE;   //!
   TBranch        *b_B0_MuMu_mu1_dr;   //!
   TBranch        *b_B0_MuMu_mu1_dr_Dimuon20_Jpsi;   //!
   TBranch        *b_B0_MuMu_mu1_dr_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_MuMu_mu1_dxysign;   //!
   TBranch        *b_B0_MuMu_mu1_dzsign;   //!
   TBranch        *b_B0_MuMu_mu2_dr;   //!
   TBranch        *b_B0_MuMu_mu2_dr_Dimuon20_Jpsi;   //!
   TBranch        *b_B0_MuMu_mu2_dr_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_MuMu_mu2_dxysign;   //!
   TBranch        *b_B0_MuMu_mu2_dzsign;   //!
   TBranch        *b_B0_MuMu_prefit_mu1_eta;   //!
   TBranch        *b_B0_MuMu_prefit_mu1_phi;   //!
   TBranch        *b_B0_MuMu_prefit_mu1_pt;   //!
   TBranch        *b_B0_MuMu_prefit_mu2_eta;   //!
   TBranch        *b_B0_MuMu_prefit_mu2_phi;   //!
   TBranch        *b_B0_MuMu_prefit_mu2_pt;   //!
   TBranch        *b_B0_MuMu_sv_prob;   //!
   TBranch        *b_B0_PVEx;   //!
   TBranch        *b_B0_PVEy;   //!
   TBranch        *b_B0_PVEz;   //!
   TBranch        *b_B0_PVx;   //!
   TBranch        *b_B0_PVy;   //!
   TBranch        *b_B0_PVz;   //!
   TBranch        *b_B0_PiPi_p1_dr_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_PiPi_p2_dr_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_PiPi_pi1_d0sig;   //!
   TBranch        *b_B0_PiPi_pi1_dxysign;   //!
   TBranch        *b_B0_PiPi_pi1_dzsign;   //!
   TBranch        *b_B0_PiPi_pi1_maxd0PV;   //!
   TBranch        *b_B0_PiPi_pi1_mind0PV;   //!
   TBranch        *b_B0_PiPi_pi2_d0sig;   //!
   TBranch        *b_B0_PiPi_pi2_dxysign;   //!
   TBranch        *b_B0_PiPi_pi2_dzsign;   //!
   TBranch        *b_B0_PiPi_pi2_maxd0PV;   //!
   TBranch        *b_B0_PiPi_pi2_mind0PV;   //!
   TBranch        *b_B0_PiPi_prefit_pi1_eta;   //!
   TBranch        *b_B0_PiPi_prefit_pi1_phi;   //!
   TBranch        *b_B0_PiPi_prefit_pi1_pt;   //!
   TBranch        *b_B0_PiPi_prefit_pi1_vx;   //!
   TBranch        *b_B0_PiPi_prefit_pi1_vy;   //!
   TBranch        *b_B0_PiPi_prefit_pi1_vz;   //!
   TBranch        *b_B0_PiPi_prefit_pi2_eta;   //!
   TBranch        *b_B0_PiPi_prefit_pi2_phi;   //!
   TBranch        *b_B0_PiPi_prefit_pi2_pt;   //!
   TBranch        *b_B0_PiPi_prefit_pi2_vx;   //!
   TBranch        *b_B0_PiPi_prefit_pi2_vy;   //!
   TBranch        *b_B0_PiPi_prefit_pi2_vz;   //!
   TBranch        *b_B0_PiPi_sv_prob;   //!
   TBranch        *b_B0_cosAlpha2D_BS;   //!
   TBranch        *b_B0_cosAlpha2D_BSwithZ;   //!
   TBranch        *b_B0_cosAlpha3D_PV;   //!
   TBranch        *b_B0_decayVtxX;   //!
   TBranch        *b_B0_decayVtxXE;   //!
   TBranch        *b_B0_decayVtxY;   //!
   TBranch        *b_B0_decayVtxYE;   //!
   TBranch        *b_B0_decayVtxZ;   //!
   TBranch        *b_B0_decayVtxZE;   //!
   TBranch        *b_B0_finalFit_JPsi_mass;   //!
   TBranch        *b_B0_finalFit_Rho_mass;   //!
   TBranch        *b_B0_finalFit_X_mass;   //!
   TBranch        *b_B0_finalFit_eta;   //!
   TBranch        *b_B0_finalFit_k0s_eta;   //!
   TBranch        *b_B0_finalFit_k0s_phi;   //!
   TBranch        *b_B0_finalFit_k0s_pt;   //!
   TBranch        *b_B0_finalFit_mass;   //!
   TBranch        *b_B0_finalFit_mu1_eta;   //!
   TBranch        *b_B0_finalFit_mu1_phi;   //!
   TBranch        *b_B0_finalFit_mu1_pt;   //!
   TBranch        *b_B0_finalFit_mu2_eta;   //!
   TBranch        *b_B0_finalFit_mu2_phi;   //!
   TBranch        *b_B0_finalFit_mu2_pt;   //!
   TBranch        *b_B0_finalFit_phi;   //!
   TBranch        *b_B0_finalFit_pi1_eta;   //!
   TBranch        *b_B0_finalFit_pi1_phi;   //!
   TBranch        *b_B0_finalFit_pi1_pt;   //!
   TBranch        *b_B0_finalFit_pi2_eta;   //!
   TBranch        *b_B0_finalFit_pi2_phi;   //!
   TBranch        *b_B0_finalFit_pi2_pt;   //!
   TBranch        *b_B0_finalFit_pt;   //!
   TBranch        *b_B0_fitted_mass_womc;   //!
   TBranch        *b_B0_lxySign_BS;   //!
   TBranch        *b_B0_lxySign_BSwithZ;   //!
   TBranch        *b_B0_lxySign_PV;   //!
   TBranch        *b_B0_svchi2;   //!
   TBranch        *b_B0_svprob;   //!
   TBranch        *b_B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_MuMu_mu1_fired_Dimuon20_Jpsi;   //!
   TBranch        *b_B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_MuMu_mu1_trackQuality;   //!
   TBranch        *b_B0_MuMu_mu2_fired_Dimuon20_Jpsi;   //!
   TBranch        *b_B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_MuMu_mu2_trackQuality;   //!
   TBranch        *b_B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_B0_dimuon_idx;   //!
   TBranch        *b_B0_dipion_idx;   //!
   TBranch        *b_B0_k0short_idx;   //!
   TBranch        *b_B0_mu1_idx;   //!
   TBranch        *b_B0_mu2_idx;   //!
   TBranch        *b_B0_pi1_idx;   //!
   TBranch        *b_B0_pi2_idx;   //!
   TBranch        *b_B0_pv3D_idx;   //!
   TBranch        *b_nJPsiToMuMu;   //!
   TBranch        *b_nK0s;   //!
   TBranch        *b_npipi;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_trackQuality;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_looseId;   //!
   TBranch        *b_Muon_softId;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_fixedGridRhoFastjetCentral;   //!
   TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_nProbeTracks;   //!
   TBranch        *b_ProbeTracks_nValidHits;   //!
   TBranch        *b_ProbeTracks_isMatchedToLooseMuon;   //!
   TBranch        *b_ProbeTracks_isMatchedToMuon;   //!
   TBranch        *b_ProbeTracks_isMatchedToSoftMuon;   //!
   TBranch        *b_HLT_Dimuon25_Jpsi;   //!
   TBranch        *b_HLT_DoubleMu4_JpsiTrk_Displaced;   //!
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_id;   //!
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_OtherPV_z;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_SV_dlen;   //!
   TBranch        *b_SV_dlenSig;   //!
   TBranch        *b_SV_dxy;   //!
   TBranch        *b_SV_dxySig;   //!
   TBranch        *b_SV_pAngle;   //!
   TBranch        *b_SV_charge;   //!
   TBranch        *b_SV_chi2;   //!
   TBranch        *b_SV_eta;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_ndof;   //!
   TBranch        *b_SV_phi;   //!
   TBranch        *b_SV_pt;   //!
   TBranch        *b_SV_x;   //!
   TBranch        *b_SV_y;   //!
   TBranch        *b_SV_z;   //!
   TBranch        *b_SV_ntracks;   //!

   B0toX3872K0s_base(TTree *tree=0, const TString channel = "X3872");
   virtual ~B0toX3872K0s_base();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   // CHANNEL
   TString channel_; // SGN = X3872 and NORM = Psi2S


   // BLINDING 
   bool  isBlind_ = true;
   float Blind_MB0_low  = 5.18;
   float Blind_MB0_high = 5.38;

   // RECONSTRUCTED PARTICLES
   int  RecoPartFillP4(const int Bidx);


   
// quadrimomenta of reco particles
   const float mMuon = 0.105658, mPion = 0.1395704, mK0s = 0.497648;
   ROOT::Math::PtEtaPhiMVector RecoP4_Mu1, RecoP4_Mu2;
   ROOT::Math::PtEtaPhiMVector RecoP4_JPsi;
   ROOT::Math::PtEtaPhiMVector RecoP4_Pi1, RecoP4_Pi2;
   ROOT::Math::PtEtaPhiMVector RecoP4_PiPi;
   ROOT::Math::PtEtaPhiMVector RecoP4_X3872;
   ROOT::Math::PtEtaPhiMVector RecoP4_K0s, RecoP4_K0sPi1, RecoP4_K0sPi2;
   ROOT::Math::PtEtaPhiMVector RecoP4_B0;

};

#endif

#ifdef B0toX3872K0s_base_cxx
B0toX3872K0s_base::B0toX3872K0s_base(TTree *tree, const TString channel) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data/xNANO_data_2023Feb06_135.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data/xNANO_data_2023Feb06_135.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);
   channel_ = channel;
}

B0toX3872K0s_base::~B0toX3872K0s_base()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t B0toX3872K0s_base::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t B0toX3872K0s_base::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void B0toX3872K0s_base::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nB0", &nB0, &b_nB0);
   fChain->SetBranchAddress("B0_BSxRaw", B0_BSxRaw, &b_B0_BSxRaw);
   fChain->SetBranchAddress("B0_BSxWithZ", B0_BSxWithZ, &b_B0_BSxWithZ);
   fChain->SetBranchAddress("B0_BSyRaw", B0_BSyRaw, &b_B0_BSyRaw);
   fChain->SetBranchAddress("B0_BSyWithZ", B0_BSyWithZ, &b_B0_BSyWithZ);
   fChain->SetBranchAddress("B0_K0_cosAlpha3D", B0_K0_cosAlpha3D, &b_B0_K0_cosAlpha3D);
   fChain->SetBranchAddress("B0_K0_lxySign_wrtBvtx", B0_K0_lxySign_wrtBvtx, &b_B0_K0_lxySign_wrtBvtx);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_D0sign", B0_K0s_matchTrack1_D0sign, &b_B0_K0s_matchTrack1_D0sign);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_dR", B0_K0s_matchTrack1_dR, &b_B0_K0s_matchTrack1_dR);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_eta", B0_K0s_matchTrack1_eta, &b_B0_K0s_matchTrack1_eta);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_maxD0Pv", B0_K0s_matchTrack1_maxD0Pv, &b_B0_K0s_matchTrack1_maxD0Pv);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_minD0Pv", B0_K0s_matchTrack1_minD0Pv, &b_B0_K0s_matchTrack1_minD0Pv);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_phi", B0_K0s_matchTrack1_phi, &b_B0_K0s_matchTrack1_phi);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_pt", B0_K0s_matchTrack1_pt, &b_B0_K0s_matchTrack1_pt);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_D0sign", B0_K0s_matchTrack2_D0sign, &b_B0_K0s_matchTrack2_D0sign);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_dR", B0_K0s_matchTrack2_dR, &b_B0_K0s_matchTrack2_dR);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_eta", B0_K0s_matchTrack2_eta, &b_B0_K0s_matchTrack2_eta);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_maxD0Pv", B0_K0s_matchTrack2_maxD0Pv, &b_B0_K0s_matchTrack2_maxD0Pv);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_minD0Pv", B0_K0s_matchTrack2_minD0Pv, &b_B0_K0s_matchTrack2_minD0Pv);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_phi", B0_K0s_matchTrack2_phi, &b_B0_K0s_matchTrack2_phi);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_pt", B0_K0s_matchTrack2_pt, &b_B0_K0s_matchTrack2_pt);
   fChain->SetBranchAddress("B0_K0s_mcFitted_eta", B0_K0s_mcFitted_eta, &b_B0_K0s_mcFitted_eta);
   fChain->SetBranchAddress("B0_K0s_mcFitted_mass", B0_K0s_mcFitted_mass, &b_B0_K0s_mcFitted_mass);
   fChain->SetBranchAddress("B0_K0s_mcFitted_phi", B0_K0s_mcFitted_phi, &b_B0_K0s_mcFitted_phi);
   fChain->SetBranchAddress("B0_K0s_mcFitted_pt", B0_K0s_mcFitted_pt, &b_B0_K0s_mcFitted_pt);
   fChain->SetBranchAddress("B0_K0s_mcFitted_svprob", B0_K0s_mcFitted_svprob, &b_B0_K0s_mcFitted_svprob);
   fChain->SetBranchAddress("B0_K0s_mcFitted_vtxX", B0_K0s_mcFitted_vtxX, &b_B0_K0s_mcFitted_vtxX);
   fChain->SetBranchAddress("B0_K0s_mcFitted_vtxXE", B0_K0s_mcFitted_vtxXE, &b_B0_K0s_mcFitted_vtxXE);
   fChain->SetBranchAddress("B0_K0s_mcFitted_vtxY", B0_K0s_mcFitted_vtxY, &b_B0_K0s_mcFitted_vtxY);
   fChain->SetBranchAddress("B0_K0s_mcFitted_vtxYE", B0_K0s_mcFitted_vtxYE, &b_B0_K0s_mcFitted_vtxYE);
   fChain->SetBranchAddress("B0_K0s_mcFitted_vtxZ", B0_K0s_mcFitted_vtxZ, &b_B0_K0s_mcFitted_vtxZ);
   fChain->SetBranchAddress("B0_K0s_mcFitted_vtxZE", B0_K0s_mcFitted_vtxZE, &b_B0_K0s_mcFitted_vtxZE);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_mass", B0_K0s_nmcFitted_mass, &b_B0_K0s_nmcFitted_mass);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_pi1eta", B0_K0s_nmcFitted_pi1eta, &b_B0_K0s_nmcFitted_pi1eta);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_pi1phi", B0_K0s_nmcFitted_pi1phi, &b_B0_K0s_nmcFitted_pi1phi);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_pi1pt", B0_K0s_nmcFitted_pi1pt, &b_B0_K0s_nmcFitted_pi1pt);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_pi2eta", B0_K0s_nmcFitted_pi2eta, &b_B0_K0s_nmcFitted_pi2eta);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_pi2phi", B0_K0s_nmcFitted_pi2phi, &b_B0_K0s_nmcFitted_pi2phi);
   fChain->SetBranchAddress("B0_K0s_nmcFitted_pi2pt", B0_K0s_nmcFitted_pi2pt, &b_B0_K0s_nmcFitted_pi2pt);
   fChain->SetBranchAddress("B0_K0s_prefit_mass", B0_K0s_prefit_mass, &b_B0_K0s_prefit_mass);
   fChain->SetBranchAddress("B0_MuMu_DCA", B0_MuMu_DCA, &b_B0_MuMu_DCA);
   fChain->SetBranchAddress("B0_MuMu_LxySign", B0_MuMu_LxySign, &b_B0_MuMu_LxySign);
   fChain->SetBranchAddress("B0_MuMu_cosAlpha", B0_MuMu_cosAlpha, &b_B0_MuMu_cosAlpha);
   fChain->SetBranchAddress("B0_MuMu_fitted_eta", B0_MuMu_fitted_eta, &b_B0_MuMu_fitted_eta);
   fChain->SetBranchAddress("B0_MuMu_fitted_mass", B0_MuMu_fitted_mass, &b_B0_MuMu_fitted_mass);
   fChain->SetBranchAddress("B0_MuMu_fitted_phi", B0_MuMu_fitted_phi, &b_B0_MuMu_fitted_phi);
   fChain->SetBranchAddress("B0_MuMu_fitted_pt", B0_MuMu_fitted_pt, &b_B0_MuMu_fitted_pt);
   fChain->SetBranchAddress("B0_MuMu_fitted_vtxX", B0_MuMu_fitted_vtxX, &b_B0_MuMu_fitted_vtxX);
   fChain->SetBranchAddress("B0_MuMu_fitted_vtxXE", B0_MuMu_fitted_vtxXE, &b_B0_MuMu_fitted_vtxXE);
   fChain->SetBranchAddress("B0_MuMu_fitted_vtxY", B0_MuMu_fitted_vtxY, &b_B0_MuMu_fitted_vtxY);
   fChain->SetBranchAddress("B0_MuMu_fitted_vtxYE", B0_MuMu_fitted_vtxYE, &b_B0_MuMu_fitted_vtxYE);
   fChain->SetBranchAddress("B0_MuMu_fitted_vtxZ", B0_MuMu_fitted_vtxZ, &b_B0_MuMu_fitted_vtxZ);
   fChain->SetBranchAddress("B0_MuMu_fitted_vtxZE", B0_MuMu_fitted_vtxZE, &b_B0_MuMu_fitted_vtxZE);
   fChain->SetBranchAddress("B0_MuMu_mu1_dr", B0_MuMu_mu1_dr, &b_B0_MuMu_mu1_dr);
   fChain->SetBranchAddress("B0_MuMu_mu1_dr_Dimuon20_Jpsi", B0_MuMu_mu1_dr_Dimuon20_Jpsi, &b_B0_MuMu_mu1_dr_Dimuon20_Jpsi);
   fChain->SetBranchAddress("B0_MuMu_mu1_dr_DoubleMu4_JpsiTrk_Displaced", B0_MuMu_mu1_dr_DoubleMu4_JpsiTrk_Displaced, &b_B0_MuMu_mu1_dr_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_MuMu_mu1_dxysign", B0_MuMu_mu1_dxysign, &b_B0_MuMu_mu1_dxysign);
   fChain->SetBranchAddress("B0_MuMu_mu1_dzsign", B0_MuMu_mu1_dzsign, &b_B0_MuMu_mu1_dzsign);
   fChain->SetBranchAddress("B0_MuMu_mu2_dr", B0_MuMu_mu2_dr, &b_B0_MuMu_mu2_dr);
   fChain->SetBranchAddress("B0_MuMu_mu2_dr_Dimuon20_Jpsi", B0_MuMu_mu2_dr_Dimuon20_Jpsi, &b_B0_MuMu_mu2_dr_Dimuon20_Jpsi);
   fChain->SetBranchAddress("B0_MuMu_mu2_dr_DoubleMu4_JpsiTrk_Displaced", B0_MuMu_mu2_dr_DoubleMu4_JpsiTrk_Displaced, &b_B0_MuMu_mu2_dr_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_MuMu_mu2_dxysign", B0_MuMu_mu2_dxysign, &b_B0_MuMu_mu2_dxysign);
   fChain->SetBranchAddress("B0_MuMu_mu2_dzsign", B0_MuMu_mu2_dzsign, &b_B0_MuMu_mu2_dzsign);
   fChain->SetBranchAddress("B0_MuMu_prefit_mu1_eta", B0_MuMu_prefit_mu1_eta, &b_B0_MuMu_prefit_mu1_eta);
   fChain->SetBranchAddress("B0_MuMu_prefit_mu1_phi", B0_MuMu_prefit_mu1_phi, &b_B0_MuMu_prefit_mu1_phi);
   fChain->SetBranchAddress("B0_MuMu_prefit_mu1_pt", B0_MuMu_prefit_mu1_pt, &b_B0_MuMu_prefit_mu1_pt);
   fChain->SetBranchAddress("B0_MuMu_prefit_mu2_eta", B0_MuMu_prefit_mu2_eta, &b_B0_MuMu_prefit_mu2_eta);
   fChain->SetBranchAddress("B0_MuMu_prefit_mu2_phi", B0_MuMu_prefit_mu2_phi, &b_B0_MuMu_prefit_mu2_phi);
   fChain->SetBranchAddress("B0_MuMu_prefit_mu2_pt", B0_MuMu_prefit_mu2_pt, &b_B0_MuMu_prefit_mu2_pt);
   fChain->SetBranchAddress("B0_MuMu_sv_prob", B0_MuMu_sv_prob, &b_B0_MuMu_sv_prob);
   fChain->SetBranchAddress("B0_PVEx", B0_PVEx, &b_B0_PVEx);
   fChain->SetBranchAddress("B0_PVEy", B0_PVEy, &b_B0_PVEy);
   fChain->SetBranchAddress("B0_PVEz", B0_PVEz, &b_B0_PVEz);
   fChain->SetBranchAddress("B0_PVx", B0_PVx, &b_B0_PVx);
   fChain->SetBranchAddress("B0_PVy", B0_PVy, &b_B0_PVy);
   fChain->SetBranchAddress("B0_PVz", B0_PVz, &b_B0_PVz);
   fChain->SetBranchAddress("B0_PiPi_p1_dr_DoubleMu4_JpsiTrk_Displaced", B0_PiPi_p1_dr_DoubleMu4_JpsiTrk_Displaced, &b_B0_PiPi_p1_dr_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_PiPi_p2_dr_DoubleMu4_JpsiTrk_Displaced", B0_PiPi_p2_dr_DoubleMu4_JpsiTrk_Displaced, &b_B0_PiPi_p2_dr_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_PiPi_pi1_d0sig", B0_PiPi_pi1_d0sig, &b_B0_PiPi_pi1_d0sig);
   fChain->SetBranchAddress("B0_PiPi_pi1_dxysign", B0_PiPi_pi1_dxysign, &b_B0_PiPi_pi1_dxysign);
   fChain->SetBranchAddress("B0_PiPi_pi1_dzsign", B0_PiPi_pi1_dzsign, &b_B0_PiPi_pi1_dzsign);
   fChain->SetBranchAddress("B0_PiPi_pi1_maxd0PV", B0_PiPi_pi1_maxd0PV, &b_B0_PiPi_pi1_maxd0PV);
   fChain->SetBranchAddress("B0_PiPi_pi1_mind0PV", B0_PiPi_pi1_mind0PV, &b_B0_PiPi_pi1_mind0PV);
   fChain->SetBranchAddress("B0_PiPi_pi2_d0sig", B0_PiPi_pi2_d0sig, &b_B0_PiPi_pi2_d0sig);
   fChain->SetBranchAddress("B0_PiPi_pi2_dxysign", B0_PiPi_pi2_dxysign, &b_B0_PiPi_pi2_dxysign);
   fChain->SetBranchAddress("B0_PiPi_pi2_dzsign", B0_PiPi_pi2_dzsign, &b_B0_PiPi_pi2_dzsign);
   fChain->SetBranchAddress("B0_PiPi_pi2_maxd0PV", B0_PiPi_pi2_maxd0PV, &b_B0_PiPi_pi2_maxd0PV);
   fChain->SetBranchAddress("B0_PiPi_pi2_mind0PV", B0_PiPi_pi2_mind0PV, &b_B0_PiPi_pi2_mind0PV);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi1_eta", B0_PiPi_prefit_pi1_eta, &b_B0_PiPi_prefit_pi1_eta);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi1_phi", B0_PiPi_prefit_pi1_phi, &b_B0_PiPi_prefit_pi1_phi);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi1_pt", B0_PiPi_prefit_pi1_pt, &b_B0_PiPi_prefit_pi1_pt);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi1_vx", B0_PiPi_prefit_pi1_vx, &b_B0_PiPi_prefit_pi1_vx);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi1_vy", B0_PiPi_prefit_pi1_vy, &b_B0_PiPi_prefit_pi1_vy);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi1_vz", B0_PiPi_prefit_pi1_vz, &b_B0_PiPi_prefit_pi1_vz);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi2_eta", B0_PiPi_prefit_pi2_eta, &b_B0_PiPi_prefit_pi2_eta);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi2_phi", B0_PiPi_prefit_pi2_phi, &b_B0_PiPi_prefit_pi2_phi);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi2_pt", B0_PiPi_prefit_pi2_pt, &b_B0_PiPi_prefit_pi2_pt);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi2_vx", B0_PiPi_prefit_pi2_vx, &b_B0_PiPi_prefit_pi2_vx);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi2_vy", B0_PiPi_prefit_pi2_vy, &b_B0_PiPi_prefit_pi2_vy);
   fChain->SetBranchAddress("B0_PiPi_prefit_pi2_vz", B0_PiPi_prefit_pi2_vz, &b_B0_PiPi_prefit_pi2_vz);
   fChain->SetBranchAddress("B0_PiPi_sv_prob", B0_PiPi_sv_prob, &b_B0_PiPi_sv_prob);
   fChain->SetBranchAddress("B0_cosAlpha2D_BS", B0_cosAlpha2D_BS, &b_B0_cosAlpha2D_BS);
   fChain->SetBranchAddress("B0_cosAlpha2D_BSwithZ", B0_cosAlpha2D_BSwithZ, &b_B0_cosAlpha2D_BSwithZ);
   fChain->SetBranchAddress("B0_cosAlpha3D_PV", B0_cosAlpha3D_PV, &b_B0_cosAlpha3D_PV);
   fChain->SetBranchAddress("B0_decayVtxX", B0_decayVtxX, &b_B0_decayVtxX);
   fChain->SetBranchAddress("B0_decayVtxXE", B0_decayVtxXE, &b_B0_decayVtxXE);
   fChain->SetBranchAddress("B0_decayVtxY", B0_decayVtxY, &b_B0_decayVtxY);
   fChain->SetBranchAddress("B0_decayVtxYE", B0_decayVtxYE, &b_B0_decayVtxYE);
   fChain->SetBranchAddress("B0_decayVtxZ", B0_decayVtxZ, &b_B0_decayVtxZ);
   fChain->SetBranchAddress("B0_decayVtxZE", B0_decayVtxZE, &b_B0_decayVtxZE);
   fChain->SetBranchAddress("B0_finalFit_JPsi_mass", B0_finalFit_JPsi_mass, &b_B0_finalFit_JPsi_mass);
   fChain->SetBranchAddress("B0_finalFit_Rho_mass", B0_finalFit_Rho_mass, &b_B0_finalFit_Rho_mass);
   fChain->SetBranchAddress("B0_finalFit_X_mass", B0_finalFit_X_mass, &b_B0_finalFit_X_mass);
   fChain->SetBranchAddress("B0_finalFit_eta", B0_finalFit_eta, &b_B0_finalFit_eta);
   fChain->SetBranchAddress("B0_finalFit_k0s_eta", B0_finalFit_k0s_eta, &b_B0_finalFit_k0s_eta);
   fChain->SetBranchAddress("B0_finalFit_k0s_phi", B0_finalFit_k0s_phi, &b_B0_finalFit_k0s_phi);
   fChain->SetBranchAddress("B0_finalFit_k0s_pt", B0_finalFit_k0s_pt, &b_B0_finalFit_k0s_pt);
   fChain->SetBranchAddress("B0_finalFit_mass", B0_finalFit_mass, &b_B0_finalFit_mass);
   fChain->SetBranchAddress("B0_finalFit_mu1_eta", B0_finalFit_mu1_eta, &b_B0_finalFit_mu1_eta);
   fChain->SetBranchAddress("B0_finalFit_mu1_phi", B0_finalFit_mu1_phi, &b_B0_finalFit_mu1_phi);
   fChain->SetBranchAddress("B0_finalFit_mu1_pt", B0_finalFit_mu1_pt, &b_B0_finalFit_mu1_pt);
   fChain->SetBranchAddress("B0_finalFit_mu2_eta", B0_finalFit_mu2_eta, &b_B0_finalFit_mu2_eta);
   fChain->SetBranchAddress("B0_finalFit_mu2_phi", B0_finalFit_mu2_phi, &b_B0_finalFit_mu2_phi);
   fChain->SetBranchAddress("B0_finalFit_mu2_pt", B0_finalFit_mu2_pt, &b_B0_finalFit_mu2_pt);
   fChain->SetBranchAddress("B0_finalFit_phi", B0_finalFit_phi, &b_B0_finalFit_phi);
   fChain->SetBranchAddress("B0_finalFit_pi1_eta", B0_finalFit_pi1_eta, &b_B0_finalFit_pi1_eta);
   fChain->SetBranchAddress("B0_finalFit_pi1_phi", B0_finalFit_pi1_phi, &b_B0_finalFit_pi1_phi);
   fChain->SetBranchAddress("B0_finalFit_pi1_pt", B0_finalFit_pi1_pt, &b_B0_finalFit_pi1_pt);
   fChain->SetBranchAddress("B0_finalFit_pi2_eta", B0_finalFit_pi2_eta, &b_B0_finalFit_pi2_eta);
   fChain->SetBranchAddress("B0_finalFit_pi2_phi", B0_finalFit_pi2_phi, &b_B0_finalFit_pi2_phi);
   fChain->SetBranchAddress("B0_finalFit_pi2_pt", B0_finalFit_pi2_pt, &b_B0_finalFit_pi2_pt);
   fChain->SetBranchAddress("B0_finalFit_pt", B0_finalFit_pt, &b_B0_finalFit_pt);
   fChain->SetBranchAddress("B0_fitted_mass_womc", B0_fitted_mass_womc, &b_B0_fitted_mass_womc);
   fChain->SetBranchAddress("B0_lxySign_BS", B0_lxySign_BS, &b_B0_lxySign_BS);
   fChain->SetBranchAddress("B0_lxySign_BSwithZ", B0_lxySign_BSwithZ, &b_B0_lxySign_BSwithZ);
   fChain->SetBranchAddress("B0_lxySign_PV", B0_lxySign_PV, &b_B0_lxySign_PV);
   fChain->SetBranchAddress("B0_svchi2", B0_svchi2, &b_B0_svchi2);
   fChain->SetBranchAddress("B0_svprob", B0_svprob, &b_B0_svprob);
   fChain->SetBranchAddress("B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced", B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced, &b_B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced", B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced, &b_B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_MuMu_mu1_fired_Dimuon20_Jpsi", B0_MuMu_mu1_fired_Dimuon20_Jpsi, &b_B0_MuMu_mu1_fired_Dimuon20_Jpsi);
   fChain->SetBranchAddress("B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced", B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced, &b_B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_MuMu_mu1_trackQuality", B0_MuMu_mu1_trackQuality, &b_B0_MuMu_mu1_trackQuality);
   fChain->SetBranchAddress("B0_MuMu_mu2_fired_Dimuon20_Jpsi", B0_MuMu_mu2_fired_Dimuon20_Jpsi, &b_B0_MuMu_mu2_fired_Dimuon20_Jpsi);
   fChain->SetBranchAddress("B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced", B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced, &b_B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_MuMu_mu2_trackQuality", B0_MuMu_mu2_trackQuality, &b_B0_MuMu_mu2_trackQuality);
   fChain->SetBranchAddress("B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced", B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced, &b_B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced", B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced, &b_B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("B0_dimuon_idx", B0_dimuon_idx, &b_B0_dimuon_idx);
   fChain->SetBranchAddress("B0_dipion_idx", B0_dipion_idx, &b_B0_dipion_idx);
   fChain->SetBranchAddress("B0_k0short_idx", B0_k0short_idx, &b_B0_k0short_idx);
   fChain->SetBranchAddress("B0_mu1_idx", B0_mu1_idx, &b_B0_mu1_idx);
   fChain->SetBranchAddress("B0_mu2_idx", B0_mu2_idx, &b_B0_mu2_idx);
   fChain->SetBranchAddress("B0_pi1_idx", B0_pi1_idx, &b_B0_pi1_idx);
   fChain->SetBranchAddress("B0_pi2_idx", B0_pi2_idx, &b_B0_pi2_idx);
   fChain->SetBranchAddress("B0_pv3D_idx", B0_pv3D_idx, &b_B0_pv3D_idx);
   fChain->SetBranchAddress("nJPsiToMuMu", &nJPsiToMuMu, &b_nJPsiToMuMu);
   fChain->SetBranchAddress("nK0s", &nK0s, &b_nK0s);
   fChain->SetBranchAddress("npipi", &npipi, &b_npipi);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_trackQuality", Muon_trackQuality, &b_Muon_trackQuality);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_looseId", Muon_looseId, &b_Muon_looseId);
   fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentral", &fixedGridRhoFastjetCentral, &b_fixedGridRhoFastjetCentral);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, &b_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("nProbeTracks", &nProbeTracks, &b_nProbeTracks);
   fChain->SetBranchAddress("ProbeTracks_nValidHits", ProbeTracks_nValidHits, &b_ProbeTracks_nValidHits);
   fChain->SetBranchAddress("ProbeTracks_isMatchedToLooseMuon", ProbeTracks_isMatchedToLooseMuon, &b_ProbeTracks_isMatchedToLooseMuon);
   fChain->SetBranchAddress("ProbeTracks_isMatchedToMuon", ProbeTracks_isMatchedToMuon, &b_ProbeTracks_isMatchedToMuon);
   fChain->SetBranchAddress("ProbeTracks_isMatchedToSoftMuon", ProbeTracks_isMatchedToSoftMuon, &b_ProbeTracks_isMatchedToSoftMuon);
   fChain->SetBranchAddress("HLT_Dimuon25_Jpsi", &HLT_Dimuon25_Jpsi, &b_HLT_Dimuon25_Jpsi);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Displaced", &HLT_DoubleMu4_JpsiTrk_Displaced, &b_HLT_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_pt", &TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", &TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", &TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_id", &TrigObj_id, &b_TrigObj_id);
   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   fChain->SetBranchAddress("SV_dxy", SV_dxy, &b_SV_dxy);
   fChain->SetBranchAddress("SV_dxySig", SV_dxySig, &b_SV_dxySig);
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   fChain->SetBranchAddress("SV_charge", SV_charge, &b_SV_charge);
   fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   fChain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   fChain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   fChain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
   fChain->SetBranchAddress("SV_ntracks", SV_ntracks, &b_SV_ntracks);
   Notify();
}

Bool_t B0toX3872K0s_base::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void B0toX3872K0s_base::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t B0toX3872K0s_base::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef B0toX3872K0s_base_cxx
