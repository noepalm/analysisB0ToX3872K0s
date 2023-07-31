void merge_trees(TString era){

    // Check if era argument is provided
    // if not, print error and exit
    if (era == "") {
        cout << "USAGE: root merge_trees.C(\"<era>\")" << endl;
        cout << "era: single, uppercase letter" << endl;
        return;
    }

    // check if era argument is a single, uppercase letter
    // if not, print error and exit
    if (era.Length() != 1 || !era.IsAlpha()) {
        cout << "ERROR: era argument must be a single, uppercase letter" << endl;
        return;
    }


    // Create TChain
    TChain *chain = new TChain("HLTemulation");

    // Open folder and iterate over all files associated to designated era
    // If match, add file to chain
    // File format is Run3_2022<era><split_number>_blind.root
    TString folder = "/eos/home-n/npalmeri/B0toX3872K0s/data/";
    TSystemDirectory dir(folder, folder);
    TList *files = dir.GetListOfFiles();

    cout << "Scanning " << folder << " for files..." << endl;

    if (files) { 
        TSystemFile *file; 
        TString fname; 
        TIter next(files); 
        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            // check if fname starts with Run3_2022 + era + any number: if so, add to chain
            if (fname.Contains(TRegexp("Run3_2022" + era + "[0-9]+_blind.root"))) {
                cout << "   + " << fname << endl;
                chain->Add(folder + fname);
            }
        } 
    }

    // If no files added to TChain, exit with error
    if (chain->GetEntries() == 0) {
        cout << "ERROR: no files found for era " << era << endl;
        return;
    }

    // If files were found, merge trees and save to file
    chain->Merge(folder + "Run3_2022" + era + "_blind.root");
    cout << "Saving merged tree as " << + "Run3_2022" + era + "_blind.root" << " to " << folder << endl;

}