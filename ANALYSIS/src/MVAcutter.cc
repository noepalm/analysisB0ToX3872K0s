#include "../include/MVAoptimizer.h"


int main (int argc, char** argv){

     // INPUT parameters
    if(argc < 4){
        std::cout << "[ERROR] USAGE : ./MVAcutter input.roo BDT-th Mpipi-th [channel] [year] "<< std::endl; 
		exit(-1);        
    }

    TString inFile(argv[1]);// = "MVA/results/UL_2017_BDTtraining_CV3.root";
    double BDT_th   = std::stod(argv[2]);
    double Mpipi_th = std::stod(argv[3]);
    TString channel("X3872");
    if(argc>4) channel = argv[4];
    int year = 2017;
    if(argc > 5) year = std::stoi(argv[5]);

    MVAoptimizer* cutter = new MVAoptimizer(inFile, BDT_th, Mpipi_th, year, channel);
    if (channel == "Psi2S"){
        cutter->set_Blind(false);
        cutter->set_LumiNorm(1.50/2.);
    }
    double Nbkg = cutter->Nbkg_extraction();
    //double Seff = cutter->Nsgn_extraction();
    cutter->makeSGNvsBKGplot();
    //cutter->Plot2D_BDToutMpipi();
    //cutter->set_BDT_cut(-20); cutter->set_Mpipi_cut(0.);
    //double Nbkg_0 = cutter->Nbkg_extraction();
    //std::cout << " SGNeff = " << Seff << std::endl;
    //std::cout << " BKG rejection " << 1- Nbkg/Nbkg_0 << std::endl;
    //cutter->makeSGNvsBKGplot();

    return 0;

}
