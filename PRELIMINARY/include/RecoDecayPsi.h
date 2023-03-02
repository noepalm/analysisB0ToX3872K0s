#ifndef RecoDecayPsi_h
#define RecoDecayPsi_h

#include "MCbase_B0toPsi2SK0s.h"

#include "TH1F.h"


class RecoDecayPsi : public MCbase_B0toPsi2SK0s {

public:
// constructor and destructor
RecoDecayPsi(TTree *tree=0, const TString & tags = "UL_MC");
virtual ~RecoDecayPsi();


// methods for the analysis
void    Loop();
void    MCtruthMatching(const bool verbose = false);
float   DeltaPT(ROOT::Math::PtEtaPhiMVector& genV, ROOT::Math::PtEtaPhiMVector& recV);
int     RecoPartFillP4(const int Bidx);
int     TriggerSelection_Muons(const int Bidx);
int     TriggerSelection_Track(const int Bidx);

private:

// [ OUTPUT ]
TString outFilePath_;
TFile*  outFile_;

// quadrimomenta of reco particles
ROOT::Math::PtEtaPhiMVector RecoP4_Mu1, RecoP4_Mu2;
ROOT::Math::PtEtaPhiMVector RecoP4_JPsi;
ROOT::Math::PtEtaPhiMVector RecoP4_Pi1, RecoP4_Pi2;
ROOT::Math::PtEtaPhiMVector RecoP4_PiPi;
ROOT::Math::PtEtaPhiMVector RecoP4_Psi2S;
ROOT::Math::PtEtaPhiMVector RecoP4_K0s, RecoP4_K0sPi1, RecoP4_K0sPi2;
ROOT::Math::PtEtaPhiMVector RecoP4_B0;

// MC truth matching ptl-idx
int MCmatch_Mum_Idx, MCmatch_Mup_Idx;
int MCmatch_Pim_Idx, MCmatch_Pip_Idx;
int MCmatch_K0s_Idx;
int MCmatch_B0_Idx;

// ... DR and DpT
float MCmatch_Mum_DRmin, MCmatch_Mup_DRmin;
float MCmatch_Mum_DpT, MCmatch_Mup_DpT;
float MCmatch_Pim_DRmin, MCmatch_Pip_DRmin;
float MCmatch_Pim_DpT, MCmatch_Pip_DpT;
float MCmatch_K0s_DRmin;
float MCmatch_K0s_DpT;
float MCmatch_B0_DRmin;
float MCmatch_B0_DpT;


}; //RecoDecayPsi

#endif
