#include "../include/B0toX3872K0s_base.h"
#include "../include/HLTapply.h"

//C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

//ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

using namespace std;


int main(int argc, char* argv[]) {

	// inputs from shell
	char inputFileName[500];
	//char outputFileName[500];
	char dataset[500];
	if ( argc < 2 ){
		std::cout << " missing argument: insert the file and the dataset you want to use :-)" << std::endl; 
		std::cout << " ./X3872Application inputFile [dataset]" << std::endl;
		return 1;
	}
	
	strcpy(inputFileName,argv[1]);
	strcpy(dataset, argv[2]);
	
	// -------------------------
	// Loading the file from a .txt
	TChain *theChain = new TChain("Events");
	char Buffer[5000];
	std::string NtupleDir;
	char MyRootFile[10000];
	TString ChainPath("");
	std::cout << " [INPUT] : " << inputFileName << std::endl;
	ifstream *inputFile = new ifstream(inputFileName);
    int Nfile = 0;

    while( !(inputFile->eof()) ){
        inputFile->getline(Buffer,500);
        if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
        {
            sscanf(Buffer,"%s",MyRootFile);
            std::cout << MyRootFile << std::endl;
            for(int i = 1; i < 100; i++){
                ChainPath = TString(MyRootFile); 
                if(ChainPath.EndsWith("_")) ChainPath.Append(Form("%d.root", i));
                else ChainPath.Append(Form("%.3d.root", i));
                
                int status = theChain->Add(TString(ChainPath));
                if (status) Nfile++; 
            }

        }
    }


	cout<<" Number of chained files : " << Nfile << std::endl; 



	inputFile->close();
	delete inputFile;

	//cout<<" Number of events: " << theChain->GetEntries()<<std::endl; 
    HLTapply HLTapp(theChain ,"prova", TString(dataset));
	HLTapp.Loop();


	return 0;
}
