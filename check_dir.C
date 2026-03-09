void check_dir() {
    TFile *f = TFile::Open("../AnalysisResultsROOTFiles/longRangeDihadronCorr/AnalysisResults_LHC25af_pass2_623297.root");
    if (!f || f->IsZombie()) {
        cout << "Cannot open file" << endl;
        return;
    }

    cout << "\n=== All directories in file ===" << endl;
    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        cout << key->GetName() << endl;
    }

    cout << "\n=== Contents of long-range-dihadron-cor_Cent_0_20_id48490 ===" << endl;
    TDirectory *dir = (TDirectory*)f->Get("long-range-dihadron-cor_Cent_0_20_id48490");
    if (dir) {
        TIter nextObj(dir->GetListOfKeys());
        TKey *objKey;
        while ((objKey = (TKey*)nextObj())) {
            cout << "  " << objKey->GetName() << " (" << objKey->GetClassName() << ")" << endl;
        }
    } else {
        cout << "Directory not found!" << endl;
    }
    
    f->Close();
}
