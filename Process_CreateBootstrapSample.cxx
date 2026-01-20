/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch)  
 * @Date: 2025-05-18 17:23:23 
 * @Last Modified by: Zhiyong Lu
 * @Last Modified time: 2025-10-24 09:57:59
 */
//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TList.h"
#include "TH2D.h"
#include "TRandom3.h"
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

    InputUnit(std::string _fileNameSuffix, Int_t _corrType, Bool_t _isNch, Bool_t _isPtdiff, Int_t _minRange, Int_t _maxRange) :
        fileNameSuffix(_fileNameSuffix), corrType(_corrType), isNch(_isNch), isPtdiff(_isPtdiff), minRange(_minRange), maxRange(_maxRange) {}
};

void CreateBootstrapSample(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange);
void CreateBootstrapSample_PtDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t pTMin, Double_t pTMax);

void Process_CreateBootstrapSample() {

    std::vector<InputUnit> inputList;
    // Ne-Ne
    // ft0a-ft0c
    // inputList.push_back(InputUnit("LHC25af_pass2_563812", kFT0AFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_554625", kFT0AFT0C, kCent, kPtDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25af_pass2_563812", kFT0AFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_554625", kFT0AFT0C, kCent, kPtDiffOff, 80, 100));
    // tpc-ft0c
    // inputList.push_back(InputUnit("LHC25af_pass2_586466", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_586466", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25af_pass2_565246", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_571364", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    //old
    inputList.push_back(InputUnit("LHC25af_pass2_565246", kTPCFT0A, kCent, kPtDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_565246", kTPCFT0A, kCent, kPtDiffOff, 80, 100));
    inputList.push_back(InputUnit("LHC25af_pass2_565246", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_565246", kTPCFT0C, kCent, kPtDiffOff, 80, 100));

    //innerring
    inputList.push_back(InputUnit("LHC25af_pass2_591387", kTPCFT0A, kCent, kPtDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_591387", kTPCFT0A, kCent, kPtDiffOff, 80, 100));
    inputList.push_back(InputUnit("LHC25af_pass2_591387", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_591387", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    //outterring
    inputList.push_back(InputUnit("LHC25af_pass2_591389", kTPCFT0A, kCent, kPtDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_591389", kTPCFT0A, kCent, kPtDiffOff, 80, 100));
    inputList.push_back(InputUnit("LHC25af_pass2_591389", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_591389", kTPCFT0C, kCent, kPtDiffOff, 80, 100));  
    // MC
    // inputList.push_back(InputUnit("LHC25af_pass2_586594", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_586594", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25af_pass2_586594", kTPCFT0A, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_586594", kTPCFT0A, kCent, kPtDiffOff, 80, 100));

    // inputList.push_back(InputUnit("LHC25af_pass2_588055", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_588055", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25af_pass2_588055", kTPCFT0A, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_588055", kTPCFT0A, kCent, kPtDiffOff, 80, 100));
    // // tpc-tpc
    // inputList.push_back(InputUnit("LHC25af_pass2_554626", kTPCTPC, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_554627", kTPCTPC, kCent, kPtDiffOff, 80, 100));
    // // inputList.push_back(InputUnit("LHC25af_pass2_554626", kTPCTPC, kCent, kPtDiffOff, 0, 20));
    // // inputList.push_back(InputUnit("LHC25af_pass2_554627", kTPCTPC, kCent, kPtDiffOff, 80, 100));
    // O-O
    // // ft0a-ft0c
    // inputList.push_back(InputUnit("LHC25ae_pass2_561907", kFT0AFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_561907", kFT0AFT0C, kCent, kPtDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_561907", kFT0AFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_561907", kFT0AFT0C, kCent, kPtDiffOff, 80, 100));
    // // tpc-ft0c
    // inputList.push_back(InputUnit("LHC25ae_pass2_587867", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_575598", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    //outterring
    // inputList.push_back(InputUnit("LHC25ae_pass2_583828", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_583826", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    //innerring
    // inputList.push_back(InputUnit("LHC25ae_pass2_583824", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_584612", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_573730", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_575598", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_572180", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_571294", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    // // tpc-ft0a
    // inputList.push_back(InputUnit("LHC25ae_pass2_587867", kTPCFT0A, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_575598", kTPCFT0A, kCent, kPtDiffOff, 80, 100));
    // //outterring
    // inputList.push_back(InputUnit("LHC25ae_pass2_583828", kTPCFT0A, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_583826", kTPCFT0A, kCent, kPtDiffOff, 80, 100));
    //innerring
    // inputList.push_back(InputUnit("LHC25ae_pass2_583824", kTPCFT0A, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_584612", kTPCFT0A, kCent, kPtDiffOff, 80, 100));

    //MC
    // inputList.push_back(InputUnit("LHC25ae_pass2_586778", kTPCFT0C, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_586778", kTPCFT0C, kCent, kPtDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_586778", kTPCFT0A, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_586778", kTPCFT0A, kCent, kPtDiffOff, 80, 100));
    
    // inputList.push_back(InputUnit("LHC25ae_pass2_573730", kTPCFT0A, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_575598", kTPCFT0A, kCent, kPtDiffOff, 80, 100));
    // inputList.push_back(InputUnit("LHC25ae_pass2_572180", kTPCFT0A, kCent, kPtDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_571294", kTPCFT0A, kCent, kPtDiffOff, 80, 100));

    for (auto input : inputList) {
        if (input.isPtdiff) {
            std::cout << "Processing Bootstrap Sample pT diff: " << input.fileNameSuffix << std::endl;
            for (int iPt = 0; iPt < pTBins.size() - 1; iPt++) {
                double pTMin = pTBins[iPt];
                double pTMax = pTBins[iPt + 1];
                CreateBootstrapSample_PtDiff(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, pTMin, pTMax);
            }
        } else {
            std::cout << "Processing Bootstrap Sample: " << input.fileNameSuffix << std::endl;
            CreateBootstrapSample(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange);
        }
    }

}

void CreateBootstrapSample(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    // 打开输入文件
    TFile* file = TFile::Open(
        Form("./ProcessOutput/Mixed_%s_%s_%i_%i_%s.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()), 
        "READ"
    );
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }

    // 读取所有样本直方图
    std::vector<TH1D*> hists;
    for (Int_t sample = 0; sample < maxSample; ++sample) {
        TH1D* h = dynamic_cast<TH1D*>(
            file->Get(Form("hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, sample))
        );
        if (!h) {
            std::cerr << "Error loading histogram for sample " << sample << std::endl;
            continue;
        }
        hists.push_back(h);
    }

    if (hists.empty()) {
        std::cerr << "No histograms loaded!" << std::endl;
        file->Close();
        return;
    }

    // 创建输出文件
    TFile* outFile = TFile::Open(
        Form("./ProcessOutput/BootstrapSample_%s_%s_%i_%i_%s.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()), 
        "RECREATE"
    );
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error creating output file!" << std::endl;
        file->Close();
        return;
    }

    TH1D* hAll = dynamic_cast<TH1D*>(file->Get(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange)));
    if (!hAll) {
        std::cerr << "Error loading all-sample histogram!" << std::endl;
        return;
    }
    hAll->SetName(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->GetXaxis()->SetTitle("#Delta#varphi");
    hAll->Write();

    // 初始化随机数生成器
    TRandom3 randGen;

    // 生成 maxSample^2 个bootstrap样本
    for (Int_t bs = 0; bs < maxSample * maxSample; ++bs) {
        // 随机选择样本索引（允许重复）
        std::vector<Int_t> selectedIndices;
        selectedIndices.clear();
        for (Int_t i = 0; i < maxSample; ++i) {
            selectedIndices.push_back(randGen.Integer(maxSample));
        }

        // 合并选中的直方图
        TH1D* hmerge = nullptr;
        Double_t totalWeight = 0.0;

        for (Int_t idx : selectedIndices) {
            if (idx < 0 || idx >= hists.size()) {
                std::cerr << "Invalid sample index: " << idx << std::endl;
                continue;
            }

            TH1D* current = hists[idx];
            Double_t weight = current->GetEntries();

            if (!hmerge) { // 第一次初始化合并直方图
                hmerge = dynamic_cast<TH1D*>(current->Clone(
                    Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs)
                ));
                hmerge->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs));
                // hmerge->Scale(weight);
                // std::cout << "bs: " << bs << " idx: " << idx << " weight: " << weight << std::endl;
                totalWeight = weight;
            } else {       // 后续累加
                // current->Scale(weight);
                // hmerge->Add(current, weight);
                hmerge->Add(current);
                totalWeight += weight;
            }
        }

        // 归一化并保存
        if (hmerge && totalWeight > 0) {
            // hmerge->Scale(1.0 / totalWeight);
            hmerge->Scale(1.0 / selectedIndices.size());
            // Printf("size: %d", selectedIndices.size());
            // hmerge->SetDirectory(outFile);
            hmerge->Write();
            // delete hmerge; // 释放内存
        }
    }

    // 清理资源
    outFile->Close();
    file->Close();
    delete outFile;
    delete file;
}

void CreateBootstrapSample_PtDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t pTMin, Double_t pTMax) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    // 打开输入文件
    TFile* file = TFile::Open(
        Form("./ProcessOutput/PtDiff/Mixed_%s_%s_%i_%i_Pt_%0.1f_%0.1f_%s.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, pTMin, pTMax, DihadronCorrTypeName[corrType].c_str()), 
        "READ"
    );
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return;
    }

    // 读取所有样本直方图
    std::vector<TH1D*> hists;
    for (Int_t sample = 0; sample < maxSample; ++sample) {
        TH1D* h = dynamic_cast<TH1D*>(
            file->Get(Form("hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, sample))
        );
        if (!h) {
            std::cerr << "Error loading histogram for sample " << sample << std::endl;
            continue;
        }
        hists.push_back(h);
    }

    if (hists.empty()) {
        std::cerr << "No histograms loaded!" << std::endl;
        file->Close();
        return;
    }

    // 创建输出文件
    TFile* outFile = TFile::Open(
        Form("./ProcessOutput/PtDiff/BootstrapSample_%s_%s_%i_%i_Pt_%0.1f_%0.1f_%s.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, pTMin, pTMax, DihadronCorrTypeName[corrType].c_str()), 
        "RECREATE"
    );
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error creating output file!" << std::endl;
        file->Close();
        return;
    }

    TH1D* hAll = dynamic_cast<TH1D*>(file->Get(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange)));
    if (!hAll) {
        std::cerr << "Error loading all-sample histogram!" << std::endl;
        return;
    }
    hAll->SetName(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->GetXaxis()->SetTitle("#Delta#varphi");
    hAll->Write();

    // 初始化随机数生成器
    TRandom3 randGen;

    // 生成 maxSample^2 个bootstrap样本
    for (Int_t bs = 0; bs < maxSample * maxSample; ++bs) {
        // 随机选择样本索引（允许重复）
        std::vector<Int_t> selectedIndices;
        selectedIndices.clear();
        for (Int_t i = 0; i < maxSample; ++i) {
            selectedIndices.push_back(randGen.Integer(maxSample));
        }

        // 合并选中的直方图
        TH1D* hmerge = nullptr;
        Double_t totalWeight = 0.0;

        for (Int_t idx : selectedIndices) {
            if (idx < 0 || idx >= hists.size()) {
                std::cerr << "Invalid sample index: " << idx << std::endl;
                continue;
            }

            TH1D* current = hists[idx];
            Double_t weight = current->GetEntries();

            if (!hmerge) { // 第一次初始化合并直方图
                hmerge = dynamic_cast<TH1D*>(current->Clone(
                    Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs)
                ));
                hmerge->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs));
                // hmerge->Scale(weight);
                // std::cout << "bs: " << bs << " idx: " << idx << " weight: " << weight << std::endl;
                totalWeight = weight;
            } else {       // 后续累加
                // current->Scale(weight);
                // hmerge->Add(current, weight);
                hmerge->Add(current);
                totalWeight += weight;
            }
        }

        // 归一化并保存
        if (hmerge && totalWeight > 0) {
            // hmerge->Scale(1.0 / totalWeight);
            hmerge->Scale(1.0 / selectedIndices.size());
            // Printf("size: %d", selectedIndices.size());
            // hmerge->SetDirectory(outFile);
            hmerge->Write();
            delete hmerge; // 释放内存
        }
    }

    // 清理资源
    outFile->Close();
    file->Close();
    delete outFile;
    delete file;
}
