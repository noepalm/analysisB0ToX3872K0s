#include "../include/MVAoptimizer.h"

// NOEMI: pass first era file as input.root
// e.g. ./MVAcutter /eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022D_BDTapplication_CV3_Psi2S.root -2 0.5 Psi2S 2022 > Psi2S_output.dat
//      ./MVAcutter /eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022D_BDTapplication_CV3_X3872.root -2 0.5 X3872 2022 > X3872_output.dat

int main (int argc, char** argv){

    // INPUT parameters
    if(argc < 4){
        std::cout << "[ERROR] USAGE : ./MVAcutter input.root BDT-th Mpipi-th [channel] [year] "<< std::endl; 
		exit(-1);        
    }

    TString inFile(argv[1]); // = /eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022D_BDTapplication_CV3_Psi2S.root
    double BDT_th   = std::stod(argv[2]);
    double Mpipi_th = std::stod(argv[3]);
    TString channel("X3872");
    if(argc>4) channel = argv[4];
    TString year("2017");
    if(argc > 5) year = TString(argv[5]);     

    // graphics option
    gStyle->SetLineScalePS(1.75);

    MVAoptimizer* cutter = new MVAoptimizer(inFile, BDT_th, Mpipi_th, year, channel);
    if (channel == "Psi2S"){
        cutter->set_Blind(false);
    }
    double Nbkg = cutter->Nbkg_extraction();
    double Nsgn = cutter->Nsgn_extraction();
    RooRealVar Nsgn_fit;
    if(channel == "Psi2S") Nsgn_fit = cutter->TotalFit_Psi2S();
    else if(channel == "X3872") Nsgn_fit = cutter->TotalFit_X3872();
    if(channel == "Psi2S") cutter->makeSGNvsBKGplot();
    std::cout << " [RESULTS] Nbkg = " << Nbkg  << " Nsgn = " << Nsgn_fit.getVal() << std::endl;

    /*** BRANCHING RATIO CALC ***/ 
    // luminosity (calc mean + error)
    std::map<std::string, std::vector<double>> lumis = {
        {"D", {2.152, 2.17, 2.156, 2.162, 2.155, 2.173, 2.16, 2.162}},
        {"E", {6.098, 6.102, 6.109, 6.1, 6.047, 6.091, 6.015, 6.059}},
        {"F", {18.179, 18.227, 18.268, 18.308, 18.275, 18.233, 18.289, 18.389}},
        {"G", {3.079, 3.086, 3.082, 3.107, 3.089, 3.134, 3.072, 3.133}}
    };

    double lumi = 0, lumi_error = 0;
    for(const auto& lumi_era : lumis){
        double sum = std::accumulate(lumi_era.second.begin(), lumi_era.second.end(), 0.0);
        double mean = sum / lumi_era.second.size();

        double sqsum = std::inner_product(lumi_era.second.begin(), lumi_era.second.end(), lumi_era.second.begin(), 0.0);
        double stdev = std::sqrt(sqsum / lumi_era.second.size() - mean * mean);
        double stddev_mean = stdev / std::sqrt(lumi_era.second.size());

        lumi += mean;
        lumi_error = std::sqrt(std::pow(lumi_error, 2) + std::pow(stddev_mean, 2));
    }

    if(channel == "Psi2S"){
        std::cout << "### PSI2S BRANCHING RATIO ###" << std::endl;

        double xsec = 5.5e11; //[fb]

        double BR_JpsiToMuMu = 0.05961; // ~6% from PDG
        double BR_JpsiToMuMu_error = 0.00033;
        double eff_tot = 1.217330e-04 * 4./7.; //from scripts/calculate_efficiencies.C
        double eff_ntot = 9913/(0.0239*0.0239);
        double eff_tot_err = std::sqrt(eff_tot*(1-eff_tot)/eff_ntot); //binomial error with poisson error on numbers

        double BR_Psi2S = Nsgn_fit.getVal()/(xsec*lumi*BR_JpsiToMuMu*eff_tot);  
        double BR_Psi2S_error = BR_Psi2S * std::sqrt( std::pow(Nsgn_fit.getError()/Nsgn_fit.getVal(),2) + std::pow(eff_tot_err/eff_tot,2) + std::pow(BR_JpsiToMuMu_error/BR_JpsiToMuMu, 2) + std::pow(lumi_error/lumi, 2) ); // error propagation as sum in quadrature of relative errors

        // print all parameters
        std::cout << Form("N = %.1f +- %.1f (rel err = %.2f%)", Nsgn_fit.getVal(), Nsgn_fit.getError(), Nsgn_fit.getError()/Nsgn_fit.getVal()*100) << std::endl;
        std::cout << "xsec = " << xsec << " fb" << std::endl;
        std::cout << Form("lumi = %.3f +- %.3f fb^-1 (rel err = %.2f%)", lumi, lumi_error, lumi_error/lumi*100) << std::endl;
        std::cout << Form("BR_JPsiToMuMu = (%.3f +- %.3f)% (rel err = %.2f%)", BR_JpsiToMuMu * 1e2, BR_JpsiToMuMu_error * 1e2, BR_JpsiToMuMu_error/BR_JpsiToMuMu*100) << std::endl;
        std::cout << Form("eff_tot = (%.3f +- %.3f) x 10^-4 (rel err = %.2f%)", eff_tot * 1e4, eff_tot_err * 1e4, eff_tot_err/eff_tot*100) << std::endl;

        std::cout << Form("BR(B0 -> Psi2S K0s) * BR(Psi2S -> JPsi Pi Pi) = (%.3f +- %.3f) x 10^-4 (rel. err. = %.2f%)", BR_Psi2S * 1e4, BR_Psi2S_error * 1e4, BR_Psi2S_error/BR_Psi2S * 100) << std::endl;

        // print reference value
        double pdg_BR_BtoPsi2S = 2.9e-4;
        double pdg_BR_BtoPsi2S_error = 0.25e-4;
        double pdg_BR_Psi2StoJPsiPiPi = 0.3468;
        double pdg_BR_Psi2StoJPsiPiPi_error = 0.003;

        double pdg_BR_tot = pdg_BR_BtoPsi2S * pdg_BR_Psi2StoJPsiPiPi;
        double pdg_BR_tot_error = pdg_BR_tot * std::sqrt(std::pow(pdg_BR_BtoPsi2S_error/pdg_BR_BtoPsi2S, 2) + std::pow(pdg_BR_Psi2StoJPsiPiPi_error/pdg_BR_Psi2StoJPsiPiPi, 2));
        std::cout << Form("[PDG: %.3f +- %.3f) x 10^-4]", pdg_BR_tot * 1e4, pdg_BR_tot_error * 1e4) << std::endl;


    } else {
        std::cout << "### X(3872) BRANCHING RATIO ###" << std::endl;

        double xsec = 5.426e11; //[fb]

        double BR_JpsiToMuMu = 0.05961; // ~6% from PDG
        double BR_JpsiToMuMu_error = 0.00033;
        double eff_tot = 2.329592e-04 * 4./7.; //from scripts/calculate_efficiencies.C
        double eff_ntot = 471983./(0.0334);
        double eff_tot_err = std::sqrt(eff_tot*(1-eff_tot)/eff_ntot); //binomial error with poisson error on numbers

        double BR_X3872 = Nsgn_fit.getVal()/(xsec*lumi*BR_JpsiToMuMu*eff_tot);  
        double BR_X3872_error = BR_X3872 * std::sqrt( std::pow(Nsgn_fit.getError()/Nsgn_fit.getVal(),2) + std::pow(eff_tot_err/eff_tot,2) + std::pow(BR_JpsiToMuMu_error/BR_JpsiToMuMu, 2) + std::pow(lumi_error/lumi, 2) ); // error propagation as sum in quadrature of relative errors

        // print all parameters
        std::cout << Form("N = %.1f +- %.1f (rel err = %.2f%)", Nsgn_fit.getVal(), Nsgn_fit.getError(), Nsgn_fit.getError()/Nsgn_fit.getVal()*100) << std::endl;
        std::cout << "xsec = " << xsec << " fb" << std::endl;
        std::cout << Form("lumi = %.3f +- %.3f fb^-1 (rel err = %.2f%)", lumi, lumi_error, lumi_error/lumi*100) << std::endl;
        std::cout << Form("BR_JPsiToMuMu = (%.3f +- %.3f)% (rel err = %.2f%)", BR_JpsiToMuMu * 1e2, BR_JpsiToMuMu_error * 1e2, BR_JpsiToMuMu_error/BR_JpsiToMuMu*100) << std::endl;
        std::cout << Form("eff_tot = (%.3f +- %.3f) x 10^-4 (rel err = %.2f%)", eff_tot * 1e4, eff_tot_err * 1e4, eff_tot_err/eff_tot*100) << std::endl;

        std::cout << Form("BR(B0 -> X3872 K0s) * BR(X3872 -> JPsi Pi Pi) = (%.3f +- %.3f) x 10^-4 (rel. err. = %.2f%)", BR_X3872 * 1e4, BR_X3872_error * 1e4, BR_X3872_error/BR_X3872 * 100) << std::endl;

        // print reference value
        double pdg_BR_BtoX = 1.1e-4 / 2.;
        double pdg_BR_BtoX_error = 0.4e-4 / 2.;
        double pdg_BR_XtoJPsiPiPi = 0.038;
        double pdg_BR_XtoJPsiPiPi_error = 0.012;

        double pdg_BR_tot = pdg_BR_BtoX * pdg_BR_XtoJPsiPiPi;
        double pdg_BR_tot_error = pdg_BR_tot * std::sqrt(std::pow(pdg_BR_BtoX_error/pdg_BR_BtoX, 2) + std::pow(pdg_BR_XtoJPsiPiPi_error/pdg_BR_XtoJPsiPiPi, 2));
        std::cout << Form("[PDG: %.3f +- %.3f) x 10^-4]", pdg_BR_tot * 1e4, pdg_BR_tot_error * 1e4) << std::endl;
    }

    //cutter->Plot2D_BDToutMpipi();
    //cutter->set_BDT_cut(-20); cutter->set_Mpipi_cut(0.);
    //double Nbkg_0 = cutter->Nbkg_extraction();
    //cutter->makeSGNvsBKGplot();
    // std::cout << " SGNeff = " << cutter->EFFsgn_extraction() << std::endl;
    //std::cout << " BKG rejection " << 1- Nbkg/Nbkg_0 << std::endl;
    //cutter->makeSGNvsBKGplot();

    return 0;

}
