#include "../include/MVAoptimizer.h"

int main (int argc, char** argv){

    // INPUT parameters
    if(argc < 4){
        std::cout << "[ERROR] USAGE : ./MVAcutter input.roo BDT-th Mpipi-th [channel] [year] "<< std::endl; 
		exit(-1);        
    }

    TString inFile(argv[1]);// = "/eos/user/c/cbasile/B0toX3872K0s/MVAresults/UL_2017_BDTapplication_CV3_X3872.root"
    double BDT_th   = std::stod(argv[2]);
    double Mpipi_th = std::stod(argv[3]);
    TString channel("X3872");
    if(argc>4) channel = argv[4];
    TString year("2017");
    if(argc > 5) year = TString(argv[5]); 

    

    MVAoptimizer* cutter = new MVAoptimizer(inFile, BDT_th, Mpipi_th, year, channel);
    if (channel == "Psi2S"){
        cutter->set_Blind(false);
    }
    double Nbkg = cutter->Nbkg_extraction();
    double Nsgn = cutter->Nsgn_extraction();
    if(channel == "Psi2S") cutter->TotalFit_Psi2S();
    cutter->makeSGNvsBKGplot();
    std::cout << " [RESULTS] Nbkg = " << Nbkg  << " Nsgn = " << Nsgn << std::endl;
    //cutter->Plot2D_BDToutMpipi();
    //cutter->set_BDT_cut(-20); cutter->set_Mpipi_cut(0.);
    //double Nbkg_0 = cutter->Nbkg_extraction();
    //cutter->makeSGNvsBKGplot();
    std::cout << " SGNeff = " << cutter->EFFsgn_extraction() << std::endl;
    //std::cout << " BKG rejection " << 1- Nbkg/Nbkg_0 << std::endl;
    //cutter->makeSGNvsBKGplot();

    return 0;

}
