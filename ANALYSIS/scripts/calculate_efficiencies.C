int count_distinct_events(TTree* tree, TEntryList* elist, TH1I* hout = NULL){
    TString leaf_name = TString(tree->GetName()) == "HLTemulation" ? "Event" : "event";

    std::map<int, int> past_events_counter;
    int nevents = 0;

    for(int i = 0; i < elist->GetN(); i++){
        tree->GetEntry(elist->GetEntry(i));
        int current_event = tree->GetLeaf(leaf_name)->GetValue();
        if(past_events_counter.find(current_event) == past_events_counter.end()){ // if current event not found already...
            nevents++; // increment distinct event counter
            past_events_counter[current_event] = 1; // save event idx in map
        } else {
            past_events_counter[current_event]++; // if current event already saved, increment counter
        }
    }

    //if output h passed, save frequency histogram
    if(hout){
        std::cout << "inside" << std::endl;
        for(auto const& pair : past_events_counter){
            hout->Fill(pair.second);
        }
    } 

    // Return number of distinct events
    return nevents;
}

int count_candidates(TTree* tree){
    int nevs = 0;

    for(int i = 0; i < tree->GetEntries(); i++){
        tree->GetEntry(i);
        nevs += tree->GetLeaf("nB0")->GetValue();
    }

    return nevs;
}

void calculate_efficiencies(TString what_to_count = "candidates", TString channel = "Psi2S"){

    // // --- SETTINGS --- //
    // TString channel = "Psi2S"; // "Psi2S" or "X3872"
    // TString what_to_count = "candidates"; // "candidates" or "events"
    // // read args
    // if(argc > 1) what_to_count = argv[1];
    // if(argc > 2) channel = argv[2];

    // --- INITIAL DATA --- //

    double efficiency_filter, efficiency_analyzer;
    if(channel == "Psi2S"){
        efficiency_filter = 0.0239;
        efficiency_analyzer = 9913./415437.;
    } else if(channel == "X3872") {
        efficiency_filter = 0.0334;
        efficiency_analyzer = 8817./471983.;
    }

    // --- CALCULATING N. OF EVENTS FOR EACH SELECTION --- //

    // 0) Retrieve initial MC events
    TFile* fMC_initial;
    if(channel == "Psi2S") fMC_initial = TFile::Open("root://xrootd-cms.infn.it///store/user/crovelli/Run32022__BuToPsi2SKs_Psi2SToJPsiPiPi_JPsiToMuMu_RhoToPiPi_TuneCP5_13p6TeV_pythia8-evtgen/BuToPsi2SKs_Psi2SToJPsiPiPi_JPsiToMuMu_RhoToPiPi.root");
    else fMC_initial = TFile::Open("root://xrootd-cms.infn.it//store/user/crovelli/Run32022__BuToX3872Ks_X3872ToJPsiRho_JPsiToMuMu_RhoToPiPi_TuneCP5_13TeV_pythia8-evtgen/BuToX3872Ks_X3872ToJPsiRho_JPsiToMuMu_RhoToPiPi.root");

    TTree* tMC_initial = (TTree*)fMC_initial->Get("Events");

    int nMC_initial;
    if(what_to_count == "candidates") nMC_initial = count_candidates(tMC_initial);
    else if(what_to_count == "events") nMC_initial = tMC_initial->GetEntries();
    else std::cout << "ERROR: what_to_count must be either 'candidates' or 'events'" << std::endl;

    // 1) Calculate HLT emulation efficiency
    TFile* fMC_HLTemul;
    if(channel == "Psi2S") fMC_HLTemul = TFile::Open("~/analysisB0ToX3872K0s/ANALYSIS/outRoot/RecoDecay_Psi2S_noMCmatch_Run3.root"); // NOTE: one entry per candidate here
    else fMC_HLTemul = TFile::Open("~/analysisB0ToX3872K0s/ANALYSIS/outRoot/RecoDecay_X3872_noMCmatch_Run3.root"); // NOTE: one entry per candidate here

    TTree* tMC_HLTemul = (TTree*)fMC_HLTemul->Get("HLTemulation");
    tMC_HLTemul->Draw(">>elist_hltemul", "", "entrylist");
    TEntryList *elist_hltemul = (TEntryList*)gDirectory->Get("elist_hltemul");

    int nMC_HLTemul;
    if(what_to_count == "candidates") nMC_HLTemul = tMC_HLTemul->GetEntries();
    else if(what_to_count == "events") nMC_HLTemul = count_distinct_events(tMC_HLTemul, elist_hltemul);


    // 2) Calculate SR selection efficiency
    TFile* fMC_BDTapplied;
    if(channel == "Psi2S") fMC_BDTapplied = TFile::Open("/eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022D_BDTapplication_CV3_Psi2S.root");
    else fMC_BDTapplied = TFile::Open("/eos/home-n/npalmeri/B0toX3872K0s/MVAresults/2022D_BDTapplication_CV3_X3872.root");

    TTree* tMC_BDTapplied = (TTree*)fMC_BDTapplied->Get("CVtraining_2022D");

    // Bounds (taken from ANALYZER/src/MVAoptimizer.C)
    double K0s_SR_low  = 0.464637;
    double K0s_SR_high = 0.530962;
    double PiPi_SR_low = 0.5;
    double PiPi_SR_high =  0.85;
    // double JpsiPiPi_SR_high = 4.0;
    // double MB_low = 5.0 , MB_high = 5.55;
    double MB_low = 5.21183 , MB_high = 5.34863;

    double JpsiPiPi_SR_low, JpsiPiPi_SR_high;
    if(channel == "Psi2S"){
        JpsiPiPi_SR_low = 3.65831;
        JpsiPiPi_SR_high= 3.71493;
    } else {
        JpsiPiPi_SR_low = 3.82562;
        JpsiPiPi_SR_high= 3.91839;
    }
    double BDTout_SR_low = -2;

    // Apply cut and count surviving events 
    
    // 2.1) Psi(2S), K0s signal region cuts
    TString PsiK0s_selection = Form("M_K0s > %f && M_K0s < %f && M_X3872 < %f && is_signal == 1", K0s_SR_low, K0s_SR_high, 4.0);
    tMC_BDTapplied->Draw(">>elist_psik0s", PsiK0s_selection, "entrylist");
    
    // count EVENTS in tMC_BDTapplied entries passing selection
    TEntryList *elist_psik0s = (TEntryList*)gDirectory->Get("elist_psik0s");

    int nMC_PsiK0scut;
    if(what_to_count == "candidates") nMC_PsiK0scut = elist_psik0s->GetN();
    else if(what_to_count == "events") nMC_PsiK0scut = count_distinct_events(tMC_BDTapplied, elist_psik0s);


    // 2.2) BDTout, PiPi signal region cuts
    TString BDToutPiPi_selection = Form("BDTout > %f && M_PiPi > %f && M_PiPi < %f", BDTout_SR_low, PiPi_SR_low, PiPi_SR_high);
    // TString BDToutPiPi_selection = Form("BDTout > %f && M_PiPi > %f && M_PiPi < %f && M_B0 > %f && M_B0 < %f", BDTout_SR_low, PiPi_SR_low, PiPi_SR_high, MB_low, MB_high);
    tMC_BDTapplied->Draw(">>elist_bdtoutpipi", PsiK0s_selection + "&&" + BDToutPiPi_selection, "entrylist");

    // count EVENTS in tMC_BDTapplied entries passing selection
    TEntryList *elist_bdtoutpipi = (TEntryList*)gDirectory->Get("elist_bdtoutpipi");

    // define frequency histogram to be saved
    TH1I* downstream_multiplicity = new TH1I("downstream_multiplicity", "downstream_multiplicity", 6, -0.5, 5.5);

    int nMC_BDToutPiPicut = count_distinct_events(tMC_BDTapplied, elist_bdtoutpipi, downstream_multiplicity); //computed even for candidate counting to get downstream multiplicity
    if(what_to_count == "candidates") nMC_BDToutPiPicut = elist_bdtoutpipi->GetN();

    // 3) signal region cuts
    TString MB0X_selection = Form("M_B0 > %f && M_B0 < %f && M_X3872 > %f && M_X3872 < %f", MB_low, MB_high, JpsiPiPi_SR_low, JpsiPiPi_SR_high);
    tMC_BDTapplied->Draw(">>elist_mb0x", PsiK0s_selection + "&&" + BDToutPiPi_selection + "&&" + MB0X_selection, "entrylist");
    TEntryList *elist_mb0x = (TEntryList*)gDirectory->Get("elist_mb0x");

    int nMC_MB0Xcut;
    if(what_to_count == "candidates") nMC_MB0Xcut = elist_mb0x->GetN();
    else if(what_to_count == "events") nMC_MB0Xcut = count_distinct_events(tMC_BDTapplied, elist_mb0x);

    double eff_tot = (double)nMC_MB0Xcut/nMC_initial * efficiency_filter * efficiency_analyzer;
    // double eff_tot = (double)nMC_BDToutPiPicut/nMC_initial * efficiency_filter * efficiency_analyzer;
    double eff_tot_err = std::sqrt(eff_tot*(1-eff_tot) / (nMC_initial/efficiency_filter/efficiency_analyzer));

    // --- PRINTING --- //
    std::cout << "\nSaved MC events = " << nMC_initial << "  ---HLT emulation--->  " << nMC_HLTemul << "  ---Psi(2S), K0s signal region--->  " << nMC_PsiK0scut << "  ---BDTout, PiPi signal region--->  " << nMC_BDToutPiPicut << " ---M_B0, M_X3872 signal region--->  " << nMC_MB0Xcut << std::endl;
    // std::cout << "\nSaved MC events = " << nMC_initial << "  ---HLT emulation--->  " << nMC_HLTemul << "  ---Psi(2S), K0s signal region--->  " << nMC_PsiK0scut << "  ---BDTout, PiPi signal region--->  " << nMC_BDToutPiPicut << std::endl;
    std::cout << "\nFilter efficiency: " << efficiency_filter * 100 << "%" << std::endl;
    std::cout << "Analyzer efficiency: " << efficiency_analyzer * 100  << "%" << std::endl;
    std::cout << "HLT emulation efficiency: " << (double)nMC_HLTemul/nMC_initial * 100 << "%" << std::endl;
    std::cout << "Psi(2S), K0s signal region efficiency: " << (double)nMC_PsiK0scut/nMC_HLTemul * 100 << "%" << std::endl;
    std::cout << "BDTout, PiPi signal region efficiency: " << (double)nMC_BDToutPiPicut/nMC_PsiK0scut * 100 << "%" << std::endl;
    std::cout << "M_B0, M_X3872 signal region efficiency: " << (double)nMC_MB0Xcut/nMC_BDToutPiPicut * 100 << "%" << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << Form("Total efficiency = (%.3f +- %.3f) x 10^-4 [%.1f%% relative error]", eff_tot*1e4, eff_tot_err*1e4, eff_tot_err/eff_tot*100) << std::endl;
    std::cout << std::fixed << setprecision(3) << std::endl; //reset format

    
    // --- DRAWING --- //
    Int_t Number = 2;
    Double_t Stops[5] = {0, .25, .50, .75, 1.};
    Double_t Red[5] = {0, 1., 1., 247./255., 61./255.};
    Double_t Green[5] = {162./255., 0.455, 240./255., 121./255., 1.};
    Double_t Blue[5] = {1., 0., 33./255., 125./255., 109./255.};
    TColor::CreateGradientColorTable(Number, Stops, Red, Green, Blue, 255, 0.5);

    std::cout << "mean multiplicity = " << downstream_multiplicity->GetMean() << " +- " << downstream_multiplicity->GetStdDev() << std::endl;
    TCanvas* c1 = new TCanvas();
    downstream_multiplicity->SetTitle("[MC] B^{0} multiplicity after selection;Multiplicity;Events");
    downstream_multiplicity->Draw("PFC PLC");
    c1->SaveAs("/eos/home-n/npalmeri/www/Analysis/Fit/MC_downstream_multiplicity.pdf");

}