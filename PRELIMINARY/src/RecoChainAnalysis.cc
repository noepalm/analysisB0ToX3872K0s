#include "../include/RecoDecayX.h"
#include "../include/RecoDecayPsi.h"

#include <iostream>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>

using namespace std;

TChain* TChainLoader(const std::string& inputFileName) {

    //================ Loading the directory path from file

    ifstream *inputFile = new ifstream(inputFileName);
    if (inputFile != nullptr) 
        std::cout << " ... [INPUT] " << inputFileName << std::endl;

    char Buffer[5000];
    char cDirPath[10000];
    int Nfiles = 0; 
    TString tree_path = "";
    const TString treeName = "/Events";

    //================ Creating chain                                                               
    TChain* chain =new TChain("Events");

    while( !(inputFile->eof()) ){
        inputFile->getline(Buffer,500);
        if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer))) sscanf(Buffer,"%s",cDirPath);
        else continue;

        std::string DirPath = std::string(cDirPath);
        std::cout << " ... Loading file from directory " << DirPath << std::endl;

        // //================ Read the directory
        struct dirent* file = NULL; 
        struct stat file_stats;
        const char* directory;
        DIR* dir_pointer = NULL;

        directory = cDirPath;
        dir_pointer = opendir(directory);//point to the directory

        while((file = readdir (dir_pointer))){
            if(file == NULL){
                std::cout << "ERROR null pointer to file" << std::endl;
                exit(-1);
            }

            if (strcmp(file->d_name, "xNANO_") < 0) continue; // skip "." and ".." and "log"

            Nfiles ++;
            //std::cout << file->d_name << std::endl;
            tree_path = DirPath + "/" + file->d_name + treeName; 
            chain->Add(tree_path);
        }

    }

    std::cout << " ... LOADING " << Nfiles << " FILES ..." << std::endl;
    cout<<" ... Total number of events: " << chain->GetEntries()<<std::endl;

    return chain;
}// TChainLoader()



int main (int argc, char* argv[]){

	if(argc < 3){
		std::cout << "... Usage ./CheckGenLevel [Indata directory] [tag] [SGN/NORM]"<< std::endl;
		return 1;
	}	

	std::string DirPath = argv[1];
	TString tag = argv[2];
	std::string channel = argv[3];
	TChain* chain = TChainLoader(DirPath);

    RecoDecayX* RecoAnalyzer = new RecoDecayX(chain, channel, tag);
	RecoAnalyzer->Loop();
	
	// if (channel == "SGN"){
	// 	RecoDecayX* RecoAnalyzer = new RecoDecayX(chain, channel, tag);
	// 	RecoAnalyzer->Loop();
	// 	delete RecoAnalyzer;
	// } else if (channel == "PROVA"){
	// 	RecoDecayPsi* RecoAnalyzer = new RecoDecayPsi(chain, tag);
	// 	RecoAnalyzer->Loop();
	// 	delete RecoAnalyzer;
	// }

	delete chain;

}//main()
