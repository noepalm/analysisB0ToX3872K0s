#ifndef HLTapply_h
#define HLTapply_h

#include "B0toX3872K0s_base.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include "TTree.h"  
#include "TCanvas.h"  
#include "TStyle.h"  
#include <iostream>
#include <fstream>
#include <stddef.h>
#include <vector>
#include <string>

#include <TMath.h>
#include "Math/Vector4D.h"
#include <Math/GenVector/VectorUtil.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include "TLorentzVector.h"



class HLTapply : public B0toX3872K0s_base{

public :

//CONSTRUCTOR - DECONSTRUCTOR
    HLTapply(TTree *tree=0, const TString outdir = "./outRoot", const TString tags = "");
    virtual ~HLTapply(){
        delete outFileTree_; 
    }

// Loop()
    void Loop();

// METTHODS
    //void GetFitParams();

    void OutTree_setup();

    int  TriggerSelection_Muons(const int Bidx);
    int  TriggerSelection_Track(const int Bidx);
    
// MEMBERS

    TString tags_;

    // I/O
    TString outFileHistoPath_;
    TFile*  outFileHisto_;
    TString outFileTreePath_;
    TFile*  outFileTree_;
    TTree*  outTree_;

    // branches (aggiungere anche pT ed eta)
    float Event, LumiBlock, Run;
    float M_B0, M_PiPi, M_X3872, M_mumu, M_K0s, M_K0sPi1, M_K0sPi2;	


    float LxySignBSz_B0, CosAlpha3DBSz_B0, pTM_B0, SVprob_B0;
    float LxySignSV_K0s;
    float SVprob_PiPi, pT_PiPi; 
    float pT_Pi1, DR_B0Pi1, D0_Pi1;

};


#endif
