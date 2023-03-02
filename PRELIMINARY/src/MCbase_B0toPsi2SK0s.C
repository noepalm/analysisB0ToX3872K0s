#define MCbase_B0toPsi2SK0s_cxx
#include "../include/MCbase_B0toPsi2SK0s.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MCbase_B0toPsi2SK0s::Loop()
{
//   In a ROOT session, you can do:
//      root> .L MCbase_B0toPsi2SK0s.C
//      root> MCbase_B0toPsi2SK0s t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
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


void MCbase_B0toPsi2SK0s::GenPartFillP4(){

	Int_t MumIdx = -1, MupIdx = -1, PimIdx = -1 , PipIdx = -1, JPsiIdx = -1, Psi2SIdx = -1, K0sIdx = -1, B0Idx = -1;
   for (UInt_t g = 0; g < nGenPart; g++){
		// MUONS
		if( (abs(GenPart_pdgId[g]) ==   isMum) &&
			 (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isJPsi) && 
			 (GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]] == isPsi2S) &&
			 (abs(GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]]]) == isB0) )
			{
			JPsiIdx = GenPart_genPartIdxMother[g];
			Psi2SIdx = GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]];
			B0Idx = GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]];
			if (GenPart_pdgId[g] ==   isMum) MumIdx = g;
			if (GenPart_pdgId[g] ==  -isMum) MupIdx = g;
			}
		// PIONS
      if( (abs(GenPart_pdgId[g]) ==   isPip) && 
			 (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isPsi2S) &&
			 (abs(GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]]) == isB0) )
			{
			if (GenPart_pdgId[g] == -isPip) PimIdx = g;
			if (GenPart_pdgId[g] ==  isPip) PipIdx = g;
			}
      if( (GenPart_pdgId[g] ==  isK0s) &&  (abs(GenPart_pdgId[GenPart_genPartIdxMother[g]]) == abs(isB0))) K0sIdx = g;


   }//on generated ptl

	if(0){
		std::cout << " IDX Mu-  / Mu+  " << MumIdx  << "\t" << MupIdx << std::endl; 
		std::cout << " IDX Pi-  / Pi+  " << PimIdx  << "\t" << PipIdx << std::endl; 
		std::cout << " IDX JPsi " << JPsiIdx << std::endl; 
		std::cout << " IDX Psi2S/ K0s  " << Psi2SIdx << "\t" << K0sIdx << std::endl; 
		std::cout << " IDX B0  " << B0Idx << std::endl; 
	}
	
	if(MumIdx > 0){   
      GenP4_Mum.SetPt(GenPart_pt[MumIdx]);   GenP4_Mum.SetEta(GenPart_eta[MumIdx]);   
		GenP4_Mum.SetPhi(GenPart_phi[MumIdx]); GenP4_Mum.SetM(mMuon);
	} else std::cout << "No generated mu- is found" << std::endl;  

	if(MupIdx > 0){   
      GenP4_Mup.SetPt(GenPart_pt[MupIdx]);   GenP4_Mup.SetEta(GenPart_eta[MupIdx]);
		GenP4_Mup.SetPhi(GenPart_phi[MupIdx]); GenP4_Mup.SetM(mMuon); 
	} else std::cout << "No generated mu+ is found" << std::endl;  

	if(JPsiIdx > 0){   
      GenP4_JPsi.SetPt(GenPart_pt[JPsiIdx]);    GenP4_JPsi.SetEta(GenPart_eta[JPsiIdx]);
		GenP4_JPsi.SetPhi(GenPart_phi[JPsiIdx]);  GenP4_JPsi.SetM(GenPart_mass[JPsiIdx]); 
	} else std::cout << "No generated JPsi is found" << std::endl;  

	if(PimIdx > 0){   
      GenP4_Pim.SetPt(GenPart_pt[PimIdx]);   GenP4_Pim.SetEta(GenPart_eta[PimIdx]);
		GenP4_Pim.SetPhi(GenPart_phi[PimIdx]); GenP4_Pim.SetM(mPion);
	} else std::cout << "No generated pi- is found" << std::endl;  

	if(PipIdx > 0){   
      GenP4_Pip.SetPt(GenPart_pt[PipIdx]);   GenP4_Pip.SetEta(GenPart_eta[PipIdx]);
		GenP4_Pip.SetPhi(GenPart_phi[PipIdx]); GenP4_Pip.SetM(mPion); 
	} else std::cout << "No generated pi+ is found" << std::endl;  

	if((PipIdx > 0) && (PimIdx > 0)){   
      GenP4_PiPi = GenP4_Pip + GenP4_Pim;
	} else std::cout << "No generated PiPi is found" << std::endl;  
	
	if(Psi2SIdx > 0){   
      GenP4_Psi2S.SetPt(GenPart_pt[Psi2SIdx]);   GenP4_Psi2S.SetEta(GenPart_eta[Psi2SIdx]);
		GenP4_Psi2S.SetPhi(GenPart_phi[Psi2SIdx]); GenP4_Psi2S.SetM(GenPart_mass[Psi2SIdx]); 
	} else std::cout << "No generated X(3872) is found" << std::endl;  

	if(K0sIdx > 0){   
      GenP4_K0s.SetPt(GenPart_pt[K0sIdx]);   GenP4_K0s.SetEta(GenPart_eta[K0sIdx]);
		GenP4_K0s.SetPhi(GenPart_phi[K0sIdx]); GenP4_K0s.SetM(GenPart_mass[K0sIdx]); 
	} else std::cout << "No generated K0s is found" << std::endl;  

	if(B0Idx >= 0){   
      GenP4_B0.SetPt(GenPart_pt[B0Idx]);   GenP4_B0.SetEta(GenPart_eta[B0Idx]);
		GenP4_B0.SetPhi(GenPart_phi[B0Idx]); GenP4_B0.SetM(GenPart_mass[B0Idx]); 
	} else std::cout << "No generated B0 is found" << std::endl;  

	
}//GenPartFillP4()

int MCbase_B0toPsi2SK0s::GenB0idx(){
	
	ROOT::Math::PtEtaPhiMVector V;
	int genB0idx = -1;
	
	for (UInt_t g = 0; g < nGenPart; g++){
		if (abs(GenPart_pdgId[g]) != isB0) continue;
		V.SetPt(GenPart_pt[g]); V.SetEta(GenPart_eta[g]); V.SetPhi(GenPart_phi[g]);
		if (ROOT::Math::VectorUtil::DeltaR(GenP4_B0, V) < 0.00001) genB0idx = g;
	}

	if (0) {
		V.SetPt(GenPart_pt[genB0idx]); V.SetEta(GenPart_eta[genB0idx]); V.SetPhi(GenPart_phi[genB0idx]);
		std::cout << "GenB0 vs found one pT " << GenP4_B0.Pt() << "\t" << V.Pt() << std::endl;
	}
	return genB0idx;
}//GenB0idx()
