#ifndef CheckGenLevPsi2S_h
#define CheckGenLevPsi2S_h

#include "MCbase_B0toPsi2SK0s.h" 

#include "TH1F.h"

class CheckGenLevPsi2S : public MCbase_B0toPsi2SK0s{

	public:
	CheckGenLevPsi2S(TTree *tree=0, const TString& tags = "MC");
	virtual ~CheckGenLevPsi2S(){ }

	void Loop();
	void GenPart_FillKinHist(ROOT::Math::PtEtaPhiMVector* GenVec, TH1* h_pt, TH1* h_eta, TH1* h_mass); 

	private:

	TString outFilePath_;




}; //CheckGenLevPsi2S

#endif
