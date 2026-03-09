#include <TFile.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <iostream>
#include <string>
#include <vector>

void DrawComparisonAndSave(TFile* fileFourier, TFile* fileTemplate, const std::string& harmonic, const std::string& outputName) {
    TH1D* hFourier = (TH1D*)fileFourier->Get(Form("%s_FourierFit", harmonic.c_str()));
    if (!hFourier) hFourier = (TH1D*)fileFourier->Get(harmonic.c_str());

    TH1D* hTemplate = (TH1D*)fileTemplate->Get(Form("%s_TemplateFit", harmonic.c_str()));
    if (!hTemplate) hTemplate = (TH1D*)fileTemplate->Get(harmonic.c_str());

    if (!hFourier || !hTemplate) {
        std::cerr << "Missing histogram " << harmonic << " in one of input files." << std::endl;
        return;
    }

    const Int_t colorFourier = kRed + 1;
    const Int_t colorTemplate = kBlue + 1;

    TCanvas* canvas = new TCanvas(Form("canvas_%s_%s", harmonic.c_str(), outputName.c_str()), Form("%s Comparison", harmonic.c_str()), 800, 600);

    Double_t maxVal = 0;
    for (Int_t i = 1; i <= hFourier->GetNbinsX(); i++) {
        Double_t val = hFourier->GetBinContent(i) + hFourier->GetBinError(i);
        if (val > maxVal) maxVal = val;
    }
    for (Int_t i = 1; i <= hTemplate->GetNbinsX(); i++) {
        Double_t val = hTemplate->GetBinContent(i) + hTemplate->GetBinError(i);
        if (val > maxVal) maxVal = val;
    }

    hFourier->SetMarkerStyle(20);
    hFourier->SetMarkerColor(colorFourier);
    hFourier->SetLineColor(colorFourier);
    hFourier->GetYaxis()->SetRangeUser(0, maxVal * 1.2);
    hFourier->Draw("E");
    hFourier->GetXaxis()->SetLimits(0, hFourier->GetXaxis()->GetXmax());
    gPad->Modified();
    gPad->Update();

    hTemplate->SetMarkerStyle(24);
    hTemplate->SetMarkerColor(colorTemplate);
    hTemplate->SetLineColor(colorTemplate);
    hTemplate->Draw("E same");

    TLegend* leg = new TLegend(0.55, 0.70, 0.90, 0.90);
    leg->AddEntry(hFourier, "FourierFit", "lep");
    leg->AddEntry(hTemplate, "TemplateFit", "lep");
    leg->Draw();

    canvas->SaveAs(outputName.c_str());
    std::cout << harmonic << " plot saved to " << outputName << std::endl;
}

void Plot_Compare3times2PC() {
    struct DatasetUnit {
        std::string runTag;
        std::string outputSuffix;
    };

    std::vector<DatasetUnit> datasets = {
        {"LHC25af_pass2_623297", ""},
        {"LHC25ae_pass2_623296", "_LHC25ae_pass2_623296"}
    };

    for (const auto& dataset : datasets) {
        TFile* fileFourier = new TFile(Form("./3times2PC/Vn_%s_kFourierFit_Cent_0_20.root", dataset.runTag.c_str()), "READ");
        TFile* fileTemplate = new TFile(Form("./3times2PC/Vn_%s_kTemplateFit_Cent_0_20.root", dataset.runTag.c_str()), "READ");

        if (!fileFourier || !fileFourier->IsOpen()) {
            std::cerr << "Cannot open FourierFit file for " << dataset.runTag << std::endl;
            continue;
        }
        if (!fileTemplate || !fileTemplate->IsOpen()) {
            std::cerr << "Cannot open TemplateFit file for " << dataset.runTag << std::endl;
            fileFourier->Close();
            delete fileFourier;
            continue;
        }

        DrawComparisonAndSave(fileFourier, fileTemplate, "hV2", Form("./3times2PC/Compare_V2_FourierFit_vs_TemplateFit%s.root", dataset.outputSuffix.c_str()));
        DrawComparisonAndSave(fileFourier, fileTemplate, "hV3", Form("./3times2PC/Compare_V3_FourierFit_vs_TemplateFit%s.root", dataset.outputSuffix.c_str()));
        DrawComparisonAndSave(fileFourier, fileTemplate, "hV4", Form("./3times2PC/Compare_V4_FourierFit_vs_TemplateFit%s.root", dataset.outputSuffix.c_str()));

        fileFourier->Close();
        fileTemplate->Close();
        delete fileFourier;
        delete fileTemplate;
    }
}
