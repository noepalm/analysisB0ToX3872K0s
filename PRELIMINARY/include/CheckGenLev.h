#ifndef CheckGenLev_h
#define CheckGenLev_h

#include "MCbase_B0toX3872K0s.h"

#include "TH1F.h"

class CheckGenLev : public MCbase_B0toX3872K0s {

    public:
    CheckGenLev(TTree *tree=0, const TString & tags = "UL_MC");
    virtual ~CheckGenLev(){ }

    void Loop();
    void GenPart_FillKinHist(ROOT::Math::PtEtaPhiMVector* GenVec, TH1* h_pt, TH1* h_eta, TH1* h_mass);

    private:

    TString outFilePath_;


};//CheckGenLev


#endif
