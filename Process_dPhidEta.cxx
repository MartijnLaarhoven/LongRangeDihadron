/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch)  
 * @Date: 2025-05-15 21:14:52 
 * @Last Modified by: Zhiyong Lu
 * @Last Modified time: 2025-10-24 14:49:02
 */
//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TList.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TProfile.h>
#include <TRandom3.h>
#include "TMath.h"
#include "THnSparse.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "./include/BasicForDihadron.h"

struct InputUnit {
    std::string fileNameSuffix;
    Int_t corrType;
    Bool_t isNch;
    Bool_t isPtdiff;
    Int_t minRange;
    Int_t maxRange;
    Bool_t isMc;

    InputUnit(std::string _fileNameSuffix, Int_t _corrType, Bool_t _isNch, Bool_t _isPtdiff, Int_t _minRange, Int_t _maxRange, Bool_t _isMc=false) :
        fileNameSuffix(_fileNameSuffix), corrType(_corrType), isNch(_isNch), isPtdiff(_isPtdiff), minRange(_minRange), maxRange(_maxRange), isMc(_isMc) {}
};

void printAxesInfo(THnSparseF* sparseHist);
void Read_dPhidEta_givenRange(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Bool_t isMc);
void Read_dPhidEta_givenRange_PtDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t pTMin, Double_t pTMax, Bool_t isMc);

// global variables
std::string collisionSystemName = "peripheral PbPb";
std::string additionalSuffix = "";

void Process_dPhidEta() {
    // 不显示窗口
    gROOT->SetBatch(kTRUE);

    std::vector<InputUnit> inputList;
    additionalSuffix = "";

    collisionSystemName = "O-O";
    // inputList.push_back(InputUnit("LHC25af_pass1_532067", kTPCFT0A, kCent, kPtDiffOff, 0, 10));
    // inputList.push_back(InputUnit("LHC25af_pass1_532068", kTPCFT0A, kCent, kPtDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25af_pass1_532067", kTPCFT0C, kCent, kPtDiffOff, 0, 10));
    // inputList.push_back(InputUnit("LHC25af_pass1_532068", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25af_pass1_531635", kFT0AFT0C, kCent, kPtDiffOff, 0, 10));
    // inputList.push_back(InputUnit("LHC25af_pass1_531635", kFT0AFT0C, kCent, kPtDiffOff, 80, 100));

    // inputList.push_back(InputUnit("LHC25af_pass1_532067", kTPCFT0A, kCent, kPtDiffOn, 0, 10));
    // inputList.push_back(InputUnit("LHC25af_pass1_532068", kTPCFT0A, kCent, kPtDiffOn, 80, 100));
    // inputList.push_back(InputUnit("LHC25af_pass1_532067", kTPCFT0C, kCent, kPtDiffOn, 0, 10));
    // inputList.push_back(InputUnit("LHC25af_pass1_532068", kTPCFT0C, kCent, kPtDiffOn, 80, 100));

    // Ne-Ne
    // ft0a-ft0c
    inputList.push_back(InputUnit("LHC25af_pass2_563812", kFT0AFT0C, kCent, kPtDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_554625", kFT0AFT0C, kCent, kPtDiffOff, 80, 100));
    inputList.push_back(InputUnit("LHC25af_pass2_563812", kFT0AFT0C, kCent, kPtDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_554625", kFT0AFT0C, kCent, kPtDiffOff, 80, 100));
    // tpc-ft0c
    inputList.push_back(InputUnit("LHC25af_pass2_565246", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_571364", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    inputList.push_back(InputUnit("LHC25af_pass2_565246", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_571364", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    // tpc-ft0a
    inputList.push_back(InputUnit("LHC25af_pass2_565246", kTPCFT0A, kCent, kPtDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_571364", kTPCFT0A, kCent, kPtDiffOff, 80, 100));
    inputList.push_back(InputUnit("LHC25af_pass2_565246", kTPCFT0A, kCent, kPtDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_571364", kTPCFT0A, kCent, kPtDiffOff, 80, 100));
    // O-O
    // inputList.push_back(InputUnit("LHC25ae_pass2_561907", kFT0AFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_561907", kFT0AFT0C, kCent, kPtDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_561907", kFT0AFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_561907", kFT0AFT0C, kCent, kPtDiffOff, 80, 100));
  
    // inputList.push_back(InputUnit("LHC25af_pass1_537547", kFT0AFT0C, kCent, kPtDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass1_537548", kFT0AFT0C, kCent, kPtDiffOn, 80, 100));
    

    for (auto input : inputList) {
        if (input.isPtdiff) {
            for (int iPt = 0; iPt < pTBins.size() - 1; iPt++) {
                double pTMin = pTBins[iPt];
                double pTMax = pTBins[iPt + 1];
                Read_dPhidEta_givenRange_PtDiff(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, pTMin, pTMax, input.isMc);
            }
        }
        else {
            Read_dPhidEta_givenRange(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, input.isMc);
        }
    }
}

void printAxesInfo(THnSparseF* sparseHist) {
    // print the information of the axes of THnSparseF
    if (!sparseHist) {
        std::cerr << "Error: Null histogram pointer!" << std::endl;
        return;
    }

    const Int_t nDims = sparseHist->GetNdimensions();
    std::cout << "Sparse histogram has " << nDims << " dimensions\n";

    for (Int_t iDim = 0; iDim < nDims; ++iDim) {
        TAxis* axis = sparseHist->GetAxis(iDim);
        if (!axis) {
        std::cerr << "Error: Axis " << iDim << " not found!" << std::endl;
        continue;
        }
        if (iDim == 0) {
            axis->SetRangeUser(0,0.001);
              // 获取当前显示的区间参数
            const int firstBin = axis->GetFirst();
            const int lastBin = axis->GetLast();
            const double firstEdge = axis->GetBinLowEdge(firstBin);
            const double lastEdge = axis->GetBinUpEdge(lastBin);

            std::cout << "=== Axis 0" << " (" << axis->GetTitle() << ") ===\n"
                        << "Current displayed range: [" << firstEdge << ", " << lastEdge << "]\n"
                        << "Using bins: " << firstBin << " (" << firstEdge << ") -> "
                        << lastBin << " (" << lastEdge << ")\n\n";
        }

        std::cout << "Axis " << iDim << ": \n"
                << "  Title: " << axis->GetTitle() << "\n"
                << "  Bins:  " << axis->GetNbins() << "\n"
                << "  Range: [" << axis->GetXmin() 
                << ", " << axis->GetXmax() << "]\n"
                << "---------------------------------\n";
    }
}

void Read_dPhidEta_givenRange(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Bool_t isMc=false) {
    TFile *file = TFile::Open(Form("../AnalysisResultsROOTFiles/longRangeDihadronCorr/AnalysisResults_%s.root", fileNameSuffix.c_str()), "READ");
    if (!file || file->IsZombie()) {
        std::cout << "Error: Cannot open file " << fileNameSuffix << std::endl;
        return;
    }

    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    // check if MCTrue folder is available
    if (isMc) {
        if (!file->Get(Form("long-range-dihadron-cor_%s%s_%d_%d/MCTrue", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange))) {
            std::cerr << "Error: MCTrue folder not found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }
    else {
        if (file->Get(Form("long-range-dihadron-cor_%s%s_%d_%d/MCTrue", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange))) {
            std::cerr << "Caution! you are using Reco or Data, but MCTrue folder is found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }


    CorrelationContainer *same = (CorrelationContainer*)file->Get(Form("long-range-dihadron-cor_%s%s_%d_%d/sameEvent_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    CorrelationContainer *mixed = (CorrelationContainer*)file->Get(Form("long-range-dihadron-cor_%s%s_%d_%d/mixedEvent_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    THnSparseD *trig = nullptr;
    if (!isMc) {
        trig = (THnSparseD*)file->Get(Form("long-range-dihadron-cor_%s%s_%d_%d/Trig_hist_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    } else {
        trig = (THnSparseD*)file->Get(Form("long-range-dihadron-cor_%s%s_%d_%d/MCTrue/MCTrig_hist_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    }
    

    if (!same || !mixed || !trig) {
        std::cerr << "Error getting histograms for " << fileNameSuffix << " with " << splitName << " and range [" << minRange << ", " << maxRange << "]" << std::endl;
        file->Close();
        delete file;
        return;
    }

    Int_t step = CorrelationContainer::kCFStepReconstructed;
    if (isMc) step = CorrelationContainer::kCFStepAll;

    THnSparseF *sparSig = (THnSparseF*)same->getPairHist()->getTHn(step);
    THnSparseF *sparMix = (THnSparseF*)mixed->getPairHist()->getTHn(step);

    sparSig->SetName(Form("sameEvent_%i_%i", minRange, maxRange));
    sparMix->SetName(Form("mixedEvent_%i_%i", minRange, maxRange));
    trig->SetName(Form("Trig_hist_%i_%i", minRange, maxRange));

    // Common axis settings for all samples
    sparSig->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(minPt+0.001, maxPt-0.001);
    sparMix->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(minPt+0.001, maxPt-0.001);
    sparSig->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(DihadrondEtaRange[corrType][0], DihadrondEtaRange[corrType][1]);
    sparMix->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(DihadrondEtaRange[corrType][0], DihadrondEtaRange[corrType][1]);
    trig->GetAxis(trigAxis_pT)->SetRangeUser(minPt+0.001, maxPt-0.001);

    // Create output file
    TFile* fout = TFile::Open(Form("./ProcessOutput/Mixed_%s%s_%s_%i_%i_%s.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), DihadronCorrTypeName[corrType].c_str()), "RECREATE");

    // Process all samples: -1 (all), 0 to maxSample-1
    for (Int_t sample = -1; sample < maxSample; ++sample) {
        // Set sample range if applicable
        if (sample >= 0) {
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            trig->GetAxis(trigAxis_sample)->SetRangeUser(sample, sample+0.001);
        } else {
            // Reset to full range for inclusive sample
            double xmin = sparSig->GetAxis(corrAxis_kSample)->GetXmin();
            double xmax = sparSig->GetAxis(corrAxis_kSample)->GetXmax();
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            xmin = sparMix->GetAxis(corrAxis_kSample)->GetXmin();
            xmax = sparMix->GetAxis(corrAxis_kSample)->GetXmax();
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            xmin = trig->GetAxis(trigAxis_sample)->GetXmin();
            xmax = trig->GetAxis(trigAxis_sample)->GetXmax();
            trig->GetAxis(trigAxis_sample)->SetRangeUser(xmin, xmax);
        }

        // Initialize variables for this sample
        TH2D* hPhiEtaSMsum = nullptr;
        TH2D* hPhiEtaSsum = nullptr;
        TH2D* hPhiEtaMsum = nullptr;
        Double_t nTriggersS = 0.;
        Int_t nz = sparSig->GetAxis(corrAxis_kVz)->GetNbins();
        // Set names with sample suffix
        TString suffix = (sample == -1) ? "" : Form("_%d", sample);

        // Vertex loop
        for (Int_t iz = 1; iz <= nz; ++iz) {
             // project the vertex axis into a 1D histogram of the trigger sparse
            TH1D* hTriggersS = (TH1D*)trig->Projection(trigAxis_Vz);
            // get the number of triggers in this vertex bin
            nTriggersS += hTriggersS->Integral(iz, iz);

            // and set the vertex range for the same and mixed event sparses
            sparSig->GetAxis(corrAxis_kVz)->SetRange(iz, iz);
            sparMix->GetAxis(corrAxis_kVz)->SetRange(iz, iz);

            //project the phi and eta axis into a 2D Histogram for both same and mixed
            TH2D *hPhiEtaS = (TH2D*)sparSig->Projection(corrAxis_kdEtaTPCTPC, corrAxis_kdPhiTPCTPC);
            TH2D *hPhiEtaM = (TH2D*)sparMix->Projection(corrAxis_kdEtaTPCTPC, corrAxis_kdPhiTPCTPC);

            // Normalize mixed event
            Double_t norm = 1.;
            Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-TMath::Pi()/2 + 0.0001);
            Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*TMath::Pi()/2 - 0.0001);
            Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(MixEventNormalizationEta[corrType]);
            Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(MixEventNormalizationEta[corrType]);
            // Clamp bin indices to valid range to avoid empty integral regions when
            // the requested normalization value lies outside the axis range.
            {
                Int_t phiN = hPhiEtaM->GetXaxis()->GetNbins();
                Int_t etaN = hPhiEtaM->GetYaxis()->GetNbins();
                Int_t ob1 = binPhi1, ob2 = binPhi2, oe1 = binEta1, oe2 = binEta2;
                binPhi1 = TMath::Max(1, TMath::Min(binPhi1, phiN));
                binPhi2 = TMath::Max(1, TMath::Min(binPhi2, phiN));
                binEta1 = TMath::Max(1, TMath::Min(binEta1, etaN));
                binEta2 = TMath::Max(1, TMath::Min(binEta2, etaN));
                if (binPhi2 < binPhi1) std::swap(binPhi1, binPhi2);
                if (binEta2 < binEta1) std::swap(binEta1, binEta2);
                if (ob1 != binPhi1 || ob2 != binPhi2 || oe1 != binEta1 || oe2 != binEta2) {
                    std::cout << "[Debug] Clamped normalization bins for " << fileNameSuffix << " range [" << minRange << "," << maxRange << "] sample " << sample << " iz " << iz << std::endl;
                }
            }
            Int_t nNormBins = (binEta2 - binEta1 + 1) * (binPhi2 - binPhi1 + 1);
            if (nNormBins <= 0) {
                norm = 0.;
            } else {
                norm = hPhiEtaM->Integral(binPhi1, binPhi2, binEta1, binEta2) / nNormBins;
            }

            // Diagnostic / safety: if norm is zero or NaN, print useful info and protect
            if (!TMath::Finite(norm) || norm == 0.) {
                double totalInt = hPhiEtaM->Integral();
                std::cout << "[Warning] zero/NaN norm for file " << fileNameSuffix << " range [" << minRange << "," << maxRange << "] sample " << sample << " iz " << iz << " corrType " << DihadronCorrTypeName[corrType] << ": norm=" << norm << " totalIntegral=" << totalInt << std::endl;
                norm = 1.0; // avoid division by zero; keeps histogram untouched
            } else if (norm < 1e-12) {
                std::cout << "[Warning] very small norm=" << norm << " for " << fileNameSuffix << " setting to 1e-12" << std::endl;
                norm = 1e-12;
            }

            // small diagnostic counts of non-zero bins (helps understand empty graphs)
            {
                Int_t nx = hPhiEtaM->GetXaxis()->GetNbins();
                Int_t ny = hPhiEtaM->GetYaxis()->GetNbins();
                Int_t nonzeroM = 0, nonzeroS = 0;
                for (Int_t ix = 1; ix <= nx; ++ix) {
                    for (Int_t iy = 1; iy <= ny; ++iy) {
                        Int_t bin = hPhiEtaM->GetBin(ix, iy);
                        if (hPhiEtaM->GetBinContent(bin) != 0.) ++nonzeroM;
                        if (hPhiEtaS->GetBinContent(bin) != 0.) ++nonzeroS;
                    }
                }
                if (nonzeroM == 0) std::cout << "[Info] mixed-event has 0 non-zero bins for " << fileNameSuffix << " sample " << sample << " iz " << iz << std::endl;
                if (nonzeroS == 0) std::cout << "[Info] same-event has 0 non-zero bins for " << fileNameSuffix << " sample " << sample << " iz " << iz << std::endl;
            }

            // Diagnostic / safety: if norm is zero or NaN, print useful info and protect
            if (!TMath::Finite(norm) || norm == 0.) {
                double totalInt = hPhiEtaM->Integral();
                std::cout << "[Warning] zero/NaN norm for file " << fileNameSuffix << " range [" << minRange << "," << maxRange << "] sample " << sample << " iz " << iz << " corrType " << DihadronCorrTypeName[corrType] << ": norm=" << norm << " totalIntegral=" << totalInt << std::endl;
                norm = 1.0; // avoid division by zero; keeps histogram untouched
            } else if (norm < 1e-12) {
                std::cout << "[Warning] very small norm=" << norm << " for " << fileNameSuffix << " setting to 1e-12" << std::endl;
                norm = 1e-12;
            }

            // small diagnostic counts of non-zero bins (helps understand empty graphs)
            {
                Int_t nx = hPhiEtaM->GetXaxis()->GetNbins();
                Int_t ny = hPhiEtaM->GetYaxis()->GetNbins();
                Int_t nonzeroM = 0, nonzeroS = 0;
                for (Int_t ix = 1; ix <= nx; ++ix) {
                    for (Int_t iy = 1; iy <= ny; ++iy) {
                        Int_t bin = hPhiEtaM->GetBin(ix, iy);
                        if (hPhiEtaM->GetBinContent(bin) != 0.) ++nonzeroM;
                        if (hPhiEtaS->GetBinContent(bin) != 0.) ++nonzeroS;
                    }
                }
                if (nonzeroM == 0) std::cout << "[Info] mixed-event has 0 non-zero bins for " << fileNameSuffix << " sample " << sample << " iz " << iz << std::endl;
                if (nonzeroS == 0) std::cout << "[Info] same-event has 0 non-zero bins for " << fileNameSuffix << " sample " << sample << " iz " << iz << std::endl;
            }

            if (!hPhiEtaMsum) {
                hPhiEtaMsum = (TH2D*)hPhiEtaM->Clone(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaMsum->Add(hPhiEtaM);
            }

            hPhiEtaM->Scale(1.0 / norm);

            TH2D* hPhiEtaSM = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SM_%d_%d_%d%s", minRange, maxRange, iz, suffix.Data()));
            // OLD: hPhiEtaSM->Divide(hPhiEtaM);
            // Protect against zero bins in denominator to avoid NaN/inf from Divide.
            // Create a safe copy of hPhiEtaM where zero bins are replaced by the
            // smallest positive bin content found (or a small epsilon if none).
            {
                TH2D* hPhiEtaM_safe = (TH2D*)hPhiEtaM->Clone(Form("%s_safe", hPhiEtaM->GetName()));
                Int_t nx = hPhiEtaM_safe->GetXaxis()->GetNbins();
                Int_t ny = hPhiEtaM_safe->GetYaxis()->GetNbins();
                double minPos = 0.;
                // find smallest positive bin content
                for (Int_t ix = 1; ix <= nx; ++ix) {
                    for (Int_t iy = 1; iy <= ny; ++iy) {
                        Int_t bin = hPhiEtaM_safe->GetBin(ix, iy);
                        double v = hPhiEtaM_safe->GetBinContent(bin);
                        if (v > 0. && (minPos == 0. || v < minPos)) minPos = v;
                    }
                }
                if (minPos == 0.) minPos = 1e-9; // fallback tiny value
                // replace zeros with minPos
                for (Int_t ix = 1; ix <= nx; ++ix) {
                    for (Int_t iy = 1; iy <= ny; ++iy) {
                        Int_t bin = hPhiEtaM_safe->GetBin(ix, iy);
                        double v = hPhiEtaM_safe->GetBinContent(bin);
                        if (v == 0.) hPhiEtaM_safe->SetBinContent(bin, minPos);
                    }
                }

                hPhiEtaSM->Divide(hPhiEtaM_safe);
                // OLD: delete hPhiEtaM_safe;  // original code freed the safe clone here
                delete hPhiEtaM_safe;
            }

            if (!hPhiEtaSMsum) {
                hPhiEtaSMsum = (TH2D*)hPhiEtaSM->Clone(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaSMsum->Add(hPhiEtaSM);
            }

            if (!hPhiEtaSsum) {
                hPhiEtaSsum = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaSsum->Add(hPhiEtaS);
            }
        }

        // Normalization and final processing
        if (nTriggersS > 0) {
            hPhiEtaSMsum->Scale(1.0 / nTriggersS);
            hPhiEtaSsum->Scale(1.0 / nTriggersS);
        }
        hPhiEtaSMsum->Scale(1.0 / hPhiEtaSMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSMsum->Scale(1.0 / hPhiEtaSMsum->GetYaxis()->GetBinWidth(1));
        TH1D* hEta = hPhiEtaSMsum->ProjectionY(Form("hEta_%d_%d%s", minRange, maxRange, suffix.Data()));
        hEta->SetTitle("#Delta#eta");
        if (corrType == kFT0AFT0C) {
            // keep phi coarse rebin but avoid additional eta rebin to get finer delta-eta bins
            hPhiEtaSMsum->Rebin2D(4, 1);
        } else {
            // keep phi rebin by 2 but no eta rebin (second arg=1) for smaller delta-eta bin width
            hPhiEtaSMsum->Rebin2D(2, 1);
        }

        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaMsum->Rebin2D(1, 1);
        if (corrType == kFT0AFT0C) {
            // don't rebin eta further for FT0A-FT0C mix to keep finer delta-eta resolution
            hPhiEtaMsum->Rebin2D(4, 1);
        }

        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaSsum->Rebin2D(1, 1);
         if (corrType == kFT0AFT0C) {
            hPhiEtaSsum->Rebin2D(4, 1);
        }

        
        hPhiEtaSMsum->SetName(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaSsum->SetName(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaMsum->SetName(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data()));

        // set the ranges
        // hPhiEtaSMsum->GetYaxis()->SetRangeUser(-1.6, 1.6);
        // hPhiEtaSsum->GetYaxis()->SetRangeUser(-1.6, 1.6);
        // hPhiEtaMsum->GetYaxis()->SetRangeUser(-1.6, 1.6);

        // Set axis titles and labels for hPhiEtaSMsum
        hPhiEtaSMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSMsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaSMsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< Centrality #leq%d", minRange, maxRange));


        // Set axis titles and labels for hPhiEtaSsum
        hPhiEtaSsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaSsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaSsum->SetTitle(Form("Same Event %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaSsum->SetTitle(Form("Same Event %d< Centrality #leq%d", minRange, maxRange));

        // Set axis titles and labels for hPhiEtaMsum
        hPhiEtaMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaMsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaMsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< Centrality #leq%d", minRange, maxRange));

        // Draw histograms
        TCanvas* c1 = new TCanvas(Form("dPhidEta %s", suffix.Data()), Form("dPhidEta %s", suffix.Data()), 1200, 800);
        c1->Divide(2, 2);

    c1->cd(1);
    TH2D* hPhiEtaSMsum_draw = (TH2D*)hPhiEtaSMsum->Clone("hPhiEtaSMsum_draw");
    // hPhiEtaSMsum_draw->GetYaxis()->SetRangeUser(-1.5, 1.5);
    hPhiEtaSMsum_draw->Draw("surf1");
    TLatex latex;
    // Determine display collision system from dataset prefix when appropriate
    std::string displayCollisionName = collisionSystemName;
    // OLD: collisionSystemName was used directly in the plot label
    if (fileNameSuffix.find("LHC25af") == 0) displayCollisionName = "Ne-Ne";
    else if (fileNameSuffix.find("LHC25ae") == 0) displayCollisionName = "O-O";
    latex.SetNDC();
    latex.SetTextSize(0.04);
    // OLD: latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
    latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", displayCollisionName.c_str(), DihadronCorrTypeName[corrType].c_str()));
    latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", minPt, maxPt, minPt, maxPt));


    c1->cd(3);
    hPhiEtaSsum->Draw("surf1");
    latex.SetNDC();
    latex.SetTextSize(0.04);
    // OLD: latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
    latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", displayCollisionName.c_str(), DihadronCorrTypeName[corrType].c_str()));
    latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", minPt, maxPt, minPt, maxPt));


    c1->cd(4);
    hPhiEtaMsum->Draw("surf1");
    latex.SetNDC();
    latex.SetTextSize(0.04);
    // OLD: latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
    latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", displayCollisionName.c_str(), DihadronCorrTypeName[corrType].c_str()));
    latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", minPt, maxPt, minPt, maxPt));

        // Saving to file
        fout->cd();
        c1->Write();
        // write canvas to file
        // c1->SaveAs(Form("./ProcessOutput/Mixed_%s_%s_%i_%i.png", fileNameSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange)));
        // delete c1;

        // Write histograms
        hPhiEtaSMsum->Write();
        hPhiEtaSsum->Write();
        hPhiEtaMsum->Write();
        hEta->Write();

        // Eta Gap processing
        TH1D* hPhiSameOverMixed_pos = hPhiEtaSMsum->ProjectionX("hPhiSameOverMixed");
        hPhiSameOverMixed_pos->SetName(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiSameOverMixed_pos->SetTitle(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiSameOverMixed_pos->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiSameOverMixed_pos->Write();

        delete hPhiEtaSMsum;
        delete hPhiEtaSMsum_draw;
        delete hPhiEtaSsum;
        delete hPhiEtaMsum;
        delete hPhiSameOverMixed_pos;
        delete c1;
        
    }

    fout->Close();
    delete fout;
    file->Close();
    delete file;
    delete sparSig;
    delete sparMix;
    std::cout << "Output file: " << Form("./ProcessOutput/Mixed_%s%s_%s_%i_%i_%s.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), DihadronCorrTypeName[corrType].c_str()) << std::endl; 
    std::cout << "Processing completed for all samples." << std::endl;
}

void Read_dPhidEta_givenRange_PtDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t pTMin, Double_t pTMax, Bool_t isMc=false) {
    if (corrType == kFT0AFT0C) {
        std::cout << "FT0A-FT0C correlation do not have pT information, return..." << std::endl;
        return;
    }
    TFile *file = TFile::Open(Form("../AnalysisResultsROOTFiles/longRangeDihadronCorr/AnalysisResults_%s.root", fileNameSuffix.c_str()), "READ");
    if (!file || file->IsZombie()) {
        std::cout << "Error: Cannot open file " << fileNameSuffix << std::endl;
        return;
    }

    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    // check if MCTrue folder is available
    if (isMc) {
        if (!file->Get(Form("long-range-dihadron-cor_%s%s_%d_%d/MCTrue", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange))) {
            std::cerr << "Error: MCTrue folder not found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }
    else {
        if (file->Get(Form("long-range-dihadron-cor_%s%s_%d_%d/MCTrue", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange))) {
            std::cerr << "Caution! you are using Reco or Data, but MCTrue folder is found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }


    CorrelationContainer *same = (CorrelationContainer*)file->Get(Form("long-range-dihadron-cor_%s%s_%d_%d/sameEvent_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    CorrelationContainer *mixed = (CorrelationContainer*)file->Get(Form("long-range-dihadron-cor_%s%s_%d_%d/mixedEvent_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    THnSparseD *trig = nullptr;
    if (!isMc) {
        trig = (THnSparseD*)file->Get(Form("long-range-dihadron-cor_%s%s_%d_%d/Trig_hist_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    } else {
        trig = (THnSparseD*)file->Get(Form("long-range-dihadron-cor_%s%s_%d_%d/MCTrue/MCTrig_hist_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    }
    

    if (!same || !mixed || !trig) {
        std::cerr << "Error getting histograms for " << fileNameSuffix << " with " << splitName << " and range [" << minRange << ", " << maxRange << "]" << std::endl;
        file->Close();
        delete file;
        return;
    }

    Int_t step = CorrelationContainer::kCFStepReconstructed;
    if (isMc) step = CorrelationContainer::kCFStepAll;

    THnSparseF *sparSig = (THnSparseF*)same->getPairHist()->getTHn(step);
    THnSparseF *sparMix = (THnSparseF*)mixed->getPairHist()->getTHn(step);

    sparSig->SetName(Form("sameEvent_%i_%i", minRange, maxRange));
    sparMix->SetName(Form("mixedEvent_%i_%i", minRange, maxRange));
    trig->SetName(Form("Trig_hist_%i_%i", minRange, maxRange));

    // Common axis settings for all samples
    // trigger particle: finer pT range
    // asscoiated particle: reference pT range 0.2~3.0 GeV/c
    sparSig->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(pTMin+0.001, pTMax-0.001);
    sparMix->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(pTMin+0.001, pTMax-0.001);
    sparSig->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(DihadrondEtaRange[corrType][0], DihadrondEtaRange[corrType][1]);
    sparMix->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(DihadrondEtaRange[corrType][0], DihadrondEtaRange[corrType][1]);
    // sparSig->GetAxis(corrAxis_kPt_TPC_asso)->SetRangeUser(minPt+0.001, maxPt-0.001);
    // sparMix->GetAxis(corrAxis_kPt_TPC_asso)->SetRangeUser(minPt+0.001, maxPt-0.001);
    trig->GetAxis(trigAxis_pT)->SetRangeUser(pTMin+0.001, pTMax-0.001);

    // Create output file
    TFile* fout = TFile::Open(Form("./ProcessOutput/PtDiff/Mixed_%s%s_%s_%i_%i_Pt_%0.1f_%0.1f_%s.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), pTMin, pTMax, DihadronCorrTypeName[corrType].c_str()), "RECREATE");

    // Process all samples: -1 (all), 0 to maxSample-1
    for (Int_t sample = -1; sample < maxSample; ++sample) {
        // Set sample range if applicable
        if (sample >= 0) {
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            trig->GetAxis(trigAxis_sample)->SetRangeUser(sample, sample+0.001);
        } else {
            // Reset to full range for inclusive sample
            double xmin = sparSig->GetAxis(corrAxis_kSample)->GetXmin();
            double xmax = sparSig->GetAxis(corrAxis_kSample)->GetXmax();
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            xmin = sparMix->GetAxis(corrAxis_kSample)->GetXmin();
            xmax = sparMix->GetAxis(corrAxis_kSample)->GetXmax();
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            xmin = trig->GetAxis(trigAxis_sample)->GetXmin();
            xmax = trig->GetAxis(trigAxis_sample)->GetXmax();
            trig->GetAxis(trigAxis_sample)->SetRangeUser(xmin, xmax);
        }

        // Initialize variables for this sample
        TH2D* hPhiEtaSMsum = nullptr;
        TH2D* hPhiEtaSsum = nullptr;
        TH2D* hPhiEtaMsum = nullptr;
        Double_t nTriggersS = 0.;
        Int_t nz = sparSig->GetAxis(corrAxis_kVz)->GetNbins();
        // Set names with sample suffix
        TString suffix = (sample == -1) ? "" : Form("_%d", sample);

        // Vertex loop
        for (Int_t iz = 1; iz <= nz; ++iz) {
             // project the vertex axis into a 1D histogram of the trigger sparse
            TH1D* hTriggersS = (TH1D*)trig->Projection(trigAxis_Vz);
            // get the number of triggers in this vertex bin
            nTriggersS += hTriggersS->Integral(iz, iz);

            // and set the vertex range for the same and mixed event sparses
            sparSig->GetAxis(corrAxis_kVz)->SetRange(iz, iz);
            sparMix->GetAxis(corrAxis_kVz)->SetRange(iz, iz);

            //project the phi and eta axis into a 2D Histogram for both same and mixed
            TH2D *hPhiEtaS = (TH2D*)sparSig->Projection(corrAxis_kdEtaTPCTPC, corrAxis_kdPhiTPCTPC);
            TH2D *hPhiEtaM = (TH2D*)sparMix->Projection(corrAxis_kdEtaTPCTPC, corrAxis_kdPhiTPCTPC);

            // Normalize mixed event
            Double_t norm = 1.;
            Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-TMath::Pi()/2 + 0.0001);
            Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*TMath::Pi()/2 - 0.0001);
            Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(MixEventNormalizationEta[corrType]);
            Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(MixEventNormalizationEta[corrType]);
            // Clamp bin indices into valid range
            {
                Int_t phiN = hPhiEtaM->GetXaxis()->GetNbins();
                Int_t etaN = hPhiEtaM->GetYaxis()->GetNbins();
                Int_t ob1 = binPhi1, ob2 = binPhi2, oe1 = binEta1, oe2 = binEta2;
                binPhi1 = TMath::Max(1, TMath::Min(binPhi1, phiN));
                binPhi2 = TMath::Max(1, TMath::Min(binPhi2, phiN));
                binEta1 = TMath::Max(1, TMath::Min(binEta1, etaN));
                binEta2 = TMath::Max(1, TMath::Min(binEta2, etaN));
                if (binPhi2 < binPhi1) std::swap(binPhi1, binPhi2);
                if (binEta2 < binEta1) std::swap(binEta1, binEta2);
                if (ob1 != binPhi1 || ob2 != binPhi2 || oe1 != binEta1 || oe2 != binEta2) {
                    std::cout << "[Debug] Clamped normalization bins for " << fileNameSuffix << " range [" << minRange << "," << maxRange << "] sample " << sample << " iz " << iz << std::endl;
                }
            }
            Int_t nNormBins = (binEta2 - binEta1 + 1) * (binPhi2 - binPhi1 + 1);
            if (nNormBins <= 0) {
                norm = 0.;
            } else {
                norm = hPhiEtaM->Integral(binPhi1, binPhi2, binEta1, binEta2) / nNormBins;
            }

            if (!hPhiEtaMsum) {
                hPhiEtaMsum = (TH2D*)hPhiEtaM->Clone(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaMsum->Add(hPhiEtaM);
            }

            hPhiEtaM->Scale(1.0 / norm);

            TH2D* hPhiEtaSM = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SM_%d_%d_%d%s", minRange, maxRange, iz, suffix.Data()));
            // OLD: hPhiEtaSM->Divide(hPhiEtaM);
            // Protect against zero bins in denominator to avoid NaN/inf from Divide.
            {
                TH2D* hPhiEtaM_safe = (TH2D*)hPhiEtaM->Clone(Form("%s_safe", hPhiEtaM->GetName()));
                Int_t nx = hPhiEtaM_safe->GetXaxis()->GetNbins();
                Int_t ny = hPhiEtaM_safe->GetYaxis()->GetNbins();
                double minPos = 0.;
                for (Int_t ix = 1; ix <= nx; ++ix) {
                    for (Int_t iy = 1; iy <= ny; ++iy) {
                        Int_t bin = hPhiEtaM_safe->GetBin(ix, iy);
                        double v = hPhiEtaM_safe->GetBinContent(bin);
                        if (v > 0. && (minPos == 0. || v < minPos)) minPos = v;
                    }
                }
                if (minPos == 0.) minPos = 1e-9; // fallback tiny value
                for (Int_t ix = 1; ix <= nx; ++ix) {
                    for (Int_t iy = 1; iy <= ny; ++iy) {
                        Int_t bin = hPhiEtaM_safe->GetBin(ix, iy);
                        double v = hPhiEtaM_safe->GetBinContent(bin);
                        if (v == 0.) hPhiEtaM_safe->SetBinContent(bin, minPos);
                    }
                }

                hPhiEtaSM->Divide(hPhiEtaM_safe);
                // OLD: delete hPhiEtaM_safe;  // original code freed the safe clone here
                delete hPhiEtaM_safe;
            }

            if (!hPhiEtaSMsum) {
                hPhiEtaSMsum = (TH2D*)hPhiEtaSM->Clone(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaSMsum->Add(hPhiEtaSM);
            }

            if (!hPhiEtaSsum) {
                hPhiEtaSsum = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaSsum->Add(hPhiEtaS);
            }
        }

        // Normalization and final processing
        if (nTriggersS > 0) {
            hPhiEtaSMsum->Scale(1.0 / nTriggersS);
            hPhiEtaSsum->Scale(1.0 / nTriggersS);
        }
        hPhiEtaSMsum->Scale(1.0 / hPhiEtaSMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSMsum->Scale(1.0 / hPhiEtaSMsum->GetYaxis()->GetBinWidth(1));
        TH1D* hEta = hPhiEtaSMsum->ProjectionY(Form("hEta_%d_%d%s", minRange, maxRange, suffix.Data()));
        hEta->SetTitle("#Delta#eta");
        if (corrType == kFT0AFT0C) {
            hPhiEtaSMsum->Rebin2D(4, 2);
        } else {
            hPhiEtaSMsum->Rebin2D(2, 2);
        }

        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaMsum->Rebin2D(1, 1);

        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaSsum->Rebin2D(1, 1);

        
        hPhiEtaSMsum->SetName(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaSsum->SetName(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaMsum->SetName(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data()));

        // set the ranges
        // hPhiEtaSMsum->GetYaxis()->SetRangeUser(-1.6, 1.6);
        // hPhiEtaSsum->GetYaxis()->SetRangeUser(-1.6, 1.6);
        // hPhiEtaMsum->GetYaxis()->SetRangeUser(-1.6, 1.6);

        // Set axis titles and labels for hPhiEtaSMsum
        hPhiEtaSMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSMsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaSMsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< Centrality #leq%d", minRange, maxRange));


        // Set axis titles and labels for hPhiEtaSsum
        hPhiEtaSsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaSsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaSsum->SetTitle(Form("Same Event %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaSsum->SetTitle(Form("Same Event %d< Centrality #leq%d", minRange, maxRange));

        // Set axis titles and labels for hPhiEtaMsum
        hPhiEtaMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaMsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaMsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< Centrality #leq%d", minRange, maxRange));

        // Draw histograms
        TCanvas* c1 = new TCanvas(Form("dPhidEta %s", suffix.Data()), Form("dPhidEta %s", suffix.Data()), 1200, 800);
        c1->Divide(2, 2);

    c1->cd(1);
    TH2D* hPhiEtaSMsum_draw = (TH2D*)hPhiEtaSMsum->Clone("hPhiEtaSMsum_draw");
    // hPhiEtaSMsum_draw->GetYaxis()->SetRangeUser(-1.5, 1.5);
    hPhiEtaSMsum_draw->Draw("surf1");
    TLatex latex;
    // Determine display collision system from dataset prefix when appropriate
    std::string displayCollisionName = collisionSystemName;
    // OLD: collisionSystemName was used directly in the plot label
    if (fileNameSuffix.find("LHC25af") == 0) displayCollisionName = "Ne-Ne";
    else if (fileNameSuffix.find("LHC25ae") == 0) displayCollisionName = "O-O";
    latex.SetNDC();
    latex.SetTextSize(0.04);
    // OLD: latex.DrawLatex(0.1, 0.90, Form("ALICE %s", collisionSystemName.c_str()));
    latex.DrawLatex(0.1, 0.90, Form("ALICE %s", displayCollisionName.c_str()));
    latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", pTMin, pTMax, minPt, maxPt));


    c1->cd(3);
    hPhiEtaSsum->Draw("surf1");
    latex.SetNDC();
    latex.SetTextSize(0.04);
    // OLD: latex.DrawLatex(0.1, 0.90, Form("ALICE %s", collisionSystemName.c_str()));
    latex.DrawLatex(0.1, 0.90, Form("ALICE %s", displayCollisionName.c_str()));
    latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", pTMin, pTMax, minPt, maxPt));


    c1->cd(4);
    hPhiEtaMsum->Draw("surf1");
    latex.SetNDC();
    latex.SetTextSize(0.04);
    // OLD: latex.DrawLatex(0.1, 0.90, Form("ALICE %s", collisionSystemName.c_str()));
    latex.DrawLatex(0.1, 0.90, Form("ALICE %s", displayCollisionName.c_str()));
    latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", pTMin, pTMax, minPt, maxPt));

        // Saving to file
        fout->cd();
        c1->Write();
        // write canvas to file
        // c1->SaveAs(Form("./ProcessOutput/Mixed_%s_%s_%i_%i.png", fileNameSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange)));
        // delete c1;

        // Write histograms
        hPhiEtaSMsum->Write();
        hPhiEtaSsum->Write();
        hPhiEtaMsum->Write();
        hEta->Write();

        // Eta Gap processing
        TH1D* hPhiSameOverMixed_pos = hPhiEtaSMsum->ProjectionX("hPhiSameOverMixed");
        hPhiSameOverMixed_pos->SetName(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiSameOverMixed_pos->SetTitle(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiSameOverMixed_pos->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiSameOverMixed_pos->Write();

        delete hPhiEtaSMsum;
        delete hPhiEtaSMsum_draw;
        delete hPhiEtaSsum;
        delete hPhiEtaMsum;
        delete hPhiSameOverMixed_pos;
        delete c1;
        
    }

    // clean up
    fout->Close();
    delete fout;
    file->Close();
    delete file;
    delete sparSig;
    delete sparMix;
    // End of process
    std::cout << "Output file: " << Form("./ProcessOutput/PtDiff/Mixed_%s%s_%s_%i_%i_Pt_%0.1f_%0.1f_%s.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), pTMin, pTMax, DihadronCorrTypeName[corrType].c_str()) << std::endl;
    std::cout << "Processing completed for all samples." << std::endl;
}

