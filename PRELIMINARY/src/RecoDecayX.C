#include "../include/RecoDecayX.h"

RecoDecayX::RecoDecayX(TTree *tree, const TString & tags) : MCbase_B0toX3872K0s (tree, tags){

    RecoP4_Mu1.SetM(mMuon); RecoP4_Mu2.SetM(mMuon);
    RecoP4_Pi1.SetM(mPion); RecoP4_Pi2.SetM(mPion);
    RecoP4_K0s.SetM(mK0s);


}//RecoDecayX()

RecoDecayX::~RecoDecayX(){};

void RecoDecayX::Loop(){

    Long64_t nentries = fChain->GetEntriesFast();
    const Long64_t Nbreak = nentries +10;
    const Long64_t Nprint = (int)(nentries/20.);

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0 || (jentry + 1 ) == Nbreak) break;
        if ((jentry+1) % Nprint == 0) std::cout << "--> " << Form("%3.0f",(float)(jentry+1)/nentries* 100.) << " \%"<< std::endl;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        // ----- GENERATOR
        GenPartFillP4();
       


    }


}//Loop()