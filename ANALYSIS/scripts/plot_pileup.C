void plot_pileup(TString infile_name, TString era){

    // EXECUTE: root 'plot_pileup.C("Run3_2022E.txt", "E")'
    
    // open txt file with list of .root data files
    ifstream infile;
    infile.open(infile_name.Data());

    // create mother histogram
    TH1F* h_pileup_tot = new TH1F("pileup_tot", "pileup_tot", 100, 0, 100);

    // read .root data file name
    TString data_list_file;
    int nlist = 0;

    while(infile >> data_list_file && nlist++ < 10) {

        std::cout << "Reading list file " << data_list_file << std::endl;

        ifstream infile_list;
        infile_list.open(("../" + data_list_file).Data());
        
        TString data_file;
        int nfiles = 0;
        while(infile_list >> data_file && nfiles++ < 10){
            std::cout << "\tReading file " << data_file << std::endl;

            TFile* f = TFile::Open(data_file);
            TTree* t = (TTree*)f->Get("Events");

            TH1F* h = new TH1F("h", "h", 100, 0, 100);
            t->Draw("B0_finalFit_mu2_pt >> h", "", "goff");
            h_pileup_tot->Add(h);

            f->Close();
        }

        infile_list.close();

    }

    // save histogram to .pdf
    std::cout << "Saving histogram to pileup_"+era+".pdf" << std::endl;

    TCanvas* c = new TCanvas("c", "c", 800, 600);
    h_pileup_tot->Draw();
    c->SaveAs("p1_" + era + ".pdf");

    // Closing
    infile.close();
}
