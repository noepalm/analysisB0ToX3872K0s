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

// read data from T2 with xrootd file by file

int main(int argc, char* argv[]) {

	// inputs from shell
	char inputFileName[1000];
	char outputDir[1000];
	char dataset_tag[1000];
    //int Nfiles = 1;
	int iFile = 1;
	if ( argc < 2 ){
		std::cout << " missing argument: insert the file and the dataset you want to use :-)" << std::endl; 
		std::cout << argv[0] << " inputFile [outpudir] [dataset-tag] [iFile] " << std::endl;
		return 1;
	}
	
	strcpy(inputFileName,argv[1]);
	strcpy(outputDir,argv[2]);
    if (argc > 3) strcpy(dataset_tag, argv[3]);
    else strcpy(dataset_tag, "default"); 
    if (argc > 4) iFile = std::stoi(argv[4]);
	else iFile = 0;
	
	// -------------------------
	// Loading the file from a .txt
	TChain *theChain = new TChain("Events");
	char Buffer[5000];
	std::string NtupleDir;
	char MyRootFile[10000];
	TString RootFileName("");
	TString ChainPath("");
	std::cout << " [INPUT] : " << inputFileName << std::endl;
	ifstream *inputFile = new ifstream(inputFileName);
    int Nfile = 0;

    while( !(inputFile->eof()) ){
        inputFile->getline(Buffer,500);
		if (strstr(Buffer,">")){
			sscanf(Buffer, "%s",MyRootFile);
			RootFileName = TString(MyRootFile);
			RootFileName.ReplaceAll(">","");
			std::cout << " file name : " << RootFileName << std::endl;
		} else if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
        {
            sscanf(Buffer,"%s",MyRootFile);
            std::cout << MyRootFile << std::endl;
            int idx_dir = iFile/1000; // maximum of 1000 files in each directory
			ChainPath = TString(MyRootFile);
			ChainPath.Append(Form("%d/" + RootFileName + "%d.root", idx_dir, iFile));
			
			theChain->Add(TString(ChainPath));
			Nfile++;
			std::cout << " + chaining " << ChainPath << std::endl; 
            
        }
    }

	cout<<" Number of chained files : " << Nfile << std::endl; 

	inputFile->close();
	delete inputFile;

	//cout<<" Number of events: " << theChain->GetEntries()<<std::endl; 
    HLTapply HLTapp(theChain ,outputDir, TString(dataset_tag));
	HLTapp.Loop();


	return 0;
}
