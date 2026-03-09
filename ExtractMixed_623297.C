// Simple macro to extract Mixed event data from AnalysisResults and create Mixed files
#include "TFile.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "./include/BasicForDihadron.h"

void ExtractMixed_623297() {
    std::string dataset = "LHC25af_pass2_623297";
    
    // Open the AnalysisResults file
    TFile* inputFile = TFile::Open(Form("../AnalysisResultsROOTFiles/longRangeDihadronCorr/AnalysisResults_%s.root", dataset.c_str()), "READ");
    
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open AnalysisResults file!" << std::endl;
        return;
    }
    
    std::cout << "Successfully opened: " << inputFile->GetName() << std::endl;
    
    // List the contents to understand structure
    std::cout << "\nFile contents:" << std::endl;
    inputFile->ls();
    
    // Close for now
    inputFile->Close();
    delete inputFile;
}
