#define MCbase_B0toX3872K0s_cxx
#include "../include/MCbase_B0toX3872K0s.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MCbase_B0toX3872K0s::Loop()
{
//   In a ROOT session, you can do:
//      root> .L MCbase_B0toX3872K0s.C
//      root> MCbase_B0toX3872K0s t
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

void MCbase_B0toX3872K0s::GenPartFillP4(){

	Int_t MumIdx = -1, MupIdx = -1, PimIdx = -1 , PipIdx = -1, JPsiIdx = -1, X3872Idx = -1, K0sIdx = -1, RhoIdx = -1, B0Idx = -1;
   for (UInt_t g = 0; g < nGenPart; g++){
		// MUONS
		if( (abs(GenPart_pdgId[g]) ==   isMum) &&
			 (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isJPsi) && 
			 (GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]] == isX3872) &&
			 (abs(GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]]]) == isB0) )
			{
			JPsiIdx = GenPart_genPartIdxMother[g];
			X3872Idx = GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]];
			B0Idx = GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]];
			if (GenPart_pdgId[g] ==   isMum) MumIdx = g;
			if (GenPart_pdgId[g] ==  -isMum) MupIdx = g;
			}
		// PIONS
      if( (abs(GenPart_pdgId[g]) ==   isPip) && 
			 (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isRho) &&
			 (GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]] == isX3872) &&
			 (abs(GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]]]) == isB0) )
			{
			RhoIdx = GenPart_genPartIdxMother[g];
			if (GenPart_pdgId[g] == -isPip) PimIdx = g;
			if (GenPart_pdgId[g] ==  isPip) PipIdx = g;
			}
      if( (GenPart_pdgId[g] ==  isK0s) &&  (abs(GenPart_pdgId[GenPart_genPartIdxMother[g]]) == abs(isB0))) K0sIdx = g;


   }//on generated ptl

	if(0){
		std::cout << " IDX Mu-  / Mu+  " << MumIdx  << "\t" << MupIdx << std::endl; 
		std::cout << " IDX Pi-  / Pi+  " << PimIdx  << "\t" << PipIdx << std::endl; 
		std::cout << " IDX JPsi / Rho  " << JPsiIdx << "\t" << RhoIdx << std::endl; 
		std::cout << " IDX X3872/ K0s  " << X3872Idx << "\t" << K0sIdx << std::endl; 
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

	if(RhoIdx > 0){   
      GenP4_Rho.SetPt(GenPart_pt[RhoIdx]);   GenP4_Rho.SetEta(GenPart_eta[RhoIdx]);
		GenP4_Rho.SetPhi(GenPart_phi[RhoIdx]); GenP4_Rho.SetM(GenPart_mass[RhoIdx]); 
	} else std::cout << "No generated rho is found" << std::endl;  
	
	if(X3872Idx > 0){   
      GenP4_X3872.SetPt(GenPart_pt[X3872Idx]);   GenP4_X3872.SetEta(GenPart_eta[X3872Idx]);
		GenP4_X3872.SetPhi(GenPart_phi[X3872Idx]); GenP4_X3872.SetM(GenPart_mass[X3872Idx]); 
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