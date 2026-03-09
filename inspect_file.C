void inspect_file() {
    TFile *f = TFile::Open("../AnalysisResultsROOTFiles/longRangeDihadronCorr/AnalysisResults_LHC25af_pass2_623297.root");
    if (!f || f->IsZombie()) {
        std::cout << "Cannot open file" << std::endl;
        return;
    }
    
    std::cout << "=== Top level contents ===" << std::endl;
    f->ls();
    
    std::cout << "\n=== Looking for long-range-dihadron directories ===" << std::endl;
    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        TString name = key->GetName();
        if (name.Contains("long-range")) {
            std::cout << "\n>>> Found directory: " << name << std::endl;
            TDirectory *dir = (TDirectory*)f->Get(name);
            if (dir) {
                std::cout << "Contents:" << std::endl;
                dir->ls();
                
                // Check for objects in this directory
                TIter nextObj(dir->GetListOfKeys());
                TKey *objKey;
                while ((objKey = (TKey*)nextObj())) {
                    std::cout << "  - " << objKey->GetName() << " (class: " << objKey->GetClassName() << ")" << std::endl;
                }
            }
        }
    }
    f->Close();
}
