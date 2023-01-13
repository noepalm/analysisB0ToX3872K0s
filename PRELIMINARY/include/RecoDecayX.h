#ifndef RecoDecayX_h
#define RecoDecayX_h

#include "../include/MCbase_B0toX3872K0s.h"


class RecoDecayX : public MCbase_B0toX3872K0s {

public:
// constructor and destructor
RecoDecayX(TTree *tree=0, const TString & tags = "UL_MC");
virtual ~RecoDecayX();


// methods for the analysis
void Loop();
int  RecoPartFillP4(const int Bidx);

private:

// quadrimomenta of reco particles
ROOT::Math::PtEtaPhiMVector RecoP4_Mu1, RecoP4_Mu2;
ROOT::Math::PtEtaPhiMVector RecoP4_JPsi;
ROOT::Math::PtEtaPhiMVector RecoP4_Pi1, RecoP4_Pi2;
ROOT::Math::PtEtaPhiMVector RecoP4_Rho;
ROOT::Math::PtEtaPhiMVector RecoP4_X3872;
ROOT::Math::PtEtaPhiMVector RecoP4_K0s, RecoP4_K0sPi1, RecoP4_K0sPi2;
ROOT::Math::PtEtaPhiMVector RecoP4_B0;


}; //RecoDecayX

#endif