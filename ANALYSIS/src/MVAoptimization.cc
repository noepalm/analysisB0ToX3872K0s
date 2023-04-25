#include "../include/MVAoptimizer.h"


#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

int main(int argc, char** argv){

    // INPUT parameters
    if(argc < 4){
        std::cout << "[ERROR] USAGE : ./MVAOptimization input.root minBDT-th minMpipi-th [year] [step-size BDTout] [step-size Mpipi]"<< std::endl; 
		exit(-1);        
    }

    
    TString inFile(argv[1]);// = "MVA/results/UL_2017_BDTtraining_CV3.root";
    double minBDT_th   = std::stod(argv[2]);
    double minMpipi_th = std::stod(argv[3]);
    int year = 2017;
    if(argc > 5) year = std::stoi(argv[4]);
    double stepBDT = 0.25;
    if(argc > 6) stepBDT = std::stod(argv[5]);
    double stepMpipi = 0.025;
    if(argc > 7) stepMpipi = std::stod(argv[6]);
    
    // slicing in BDTout and Mpipi
    const double MaxBDTout = 2.1;
    const double MaxMpipi  = 0.9;
    const int NpointsBDT   = (int)((MaxBDTout - minBDT_th) / stepBDT) + 1;
    const int NpointsMpipi = (int)((MaxMpipi  - minMpipi_th) / stepMpipi) + 1;
    std::cout << Form(" --> starting optimization with %d / %d points along BDTout / Mpipi", NpointsBDT, NpointsMpipi) << std::endl;
    
    // TTree output
    TString OutFilePath = Form("./outRoot/PunziOptimization/UL_%d_PunziOptimization", year);
	if( ( stepBDT < 0.05 ) || ( stepMpipi < 0.01 ) ) OutFilePath.Append("Fine"); 
    else OutFilePath.Append("Rough"); 
	OutFilePath.Append(".root");

	double Bs, Ss, Psign, Psign_error;
    double BDT_th, Mpipi_th;
	TFile* OutFile = new TFile(OutFilePath, "RECREATE");
	TTree* OutTree = new TTree("CutOpt", "");
	OutTree->Branch("BDT_th", &BDT_th, "BDT_th/D");
	OutTree->Branch("Mpipi_th", &Mpipi_th, "Mpipi_th/D");
	//OutTree->Branch("nB_Sreg", &Bs, "Bs/D");
	//OutTree->Branch("nS_Sreg", &Ss, "Ss/D");
	OutTree->Branch("SignPunzi", &Psign, "SignPunzi/D");
	OutTree->Branch("SignPunzi_error", &Psign_error, "SignPunzi_error/D");
	
    MVAoptimizer* opt = new MVAoptimizer(inFile, minBDT_th, minMpipi_th, year);

	std::cout << "****  START OPTIMIZATION ****\n" << std::endl;
    
    for (int kx = 0; kx < NpointsBDT; kx++){
		BDT_th = minBDT_th + kx * stepBDT;
		if (BDT_th > MaxBDTout) BDT_th = MaxBDTout;
		for (int im = 0; im < NpointsMpipi; im++){
            Mpipi_th = minMpipi_th + im * stepMpipi;
			if(Mpipi_th > MaxMpipi) Mpipi_th = MaxMpipi;
            opt->set_BDT_cut(BDT_th); opt->set_Mpipi_cut(Mpipi_th);
            Psign = opt->PunziSign(&Psign_error);
            //Psign = 0; Psign_error = 0;
            std::cout<< Form(" = Punzi significance %f +/- %f  [BDTout > %f & Mpipi > %f]"
                        , Psign, Psign_error, BDT_th, Mpipi_th) << std::endl;
            OutTree->Fill();
        }
    }
    std::cout << "****  END OPTIMIZATION ****\n" << std::endl;

    // save aoutput
    OutFile->cd();
	OutTree->Write();
	std::cout << " \n\nWRITTEN TREE \"" << OutTree->GetName() << "\" IN FILE " << OutFilePath << std::endl;
	
	OutFile->Close();
	return 0;

}
