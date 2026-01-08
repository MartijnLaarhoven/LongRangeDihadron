/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch)  
 * @Date: 2025-05-18 13:51:01 
 * @Last Modified by: Zhiyong Lu
 * @Last Modified time: 2025-10-24 17:03:59
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
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "./include/BasicForDihadron.h"
#include "./include/Bootstrap.h"
#include "./include/ErrorPropagation.h"
#include "./include/plotting.h"

// define struct
struct InputUnit {
    std::string fileNameSuffix;
    Int_t corrType;
    Int_t method;
    Int_t minRange;
    Int_t maxRange;
    InputUnit(std::string _fileNameSuffix, Int_t _corrType, Int_t _method, Int_t _minRange, Int_t _maxRange) :
        fileNameSuffix(_fileNameSuffix), corrType(_corrType), method(_method), minRange(_minRange), maxRange(_maxRange) {}
};

struct ConfigUnit {
    Bool_t isNch;
    Bool_t isPtDiff;
    std::vector<InputUnit> dataList;
    std::string outputFileName;
    Bool_t constructed;
    ConfigUnit(Bool_t _isNch, Bool_t _isPtDiff,std::vector<InputUnit> _dataList, std::string _outputFileName) :
        isNch(_isNch), isPtDiff(_isPtDiff), dataList(_dataList), outputFileName(_outputFileName) {
            constructed = false;
            if (_dataList.size() != 3) {
                constructed = false;
                Printf("The input vector don't have a size of 3!");
                return;
            }
            // Allow arbitrary three-input combinations (e.g. TPC_EtaNeg×FT0C, TPC_EtaPos×FT0C, TPC_EtaNeg×TPC_EtaPos)
            // OLD: enforced kTPCFT0A, kTPCFT0C, kFT0AFT0C ordering and types.
            if (!(_dataList[0].minRange == _dataList[1].minRange && _dataList[0].minRange == _dataList[2].minRange && _dataList[0].maxRange == _dataList[1].maxRange && _dataList[0].maxRange == _dataList[2].maxRange)) {
                Printf("Warning: you are not using the same Cent/Nch range for one config:");
                Printf("kTPCFT0A %d %d, kTPCFT0C %d %d, kFT0AFT0C %d %d", _dataList[0].minRange, _dataList[0].maxRange, _dataList[1].minRange, _dataList[1].maxRange, _dataList[2].minRange, _dataList[2].maxRange);
            }
            constructed = true;
        }
};

struct VnUnit {
    Double_t v2;
    Double_t v2_err;
    Double_t v3;
    Double_t v3_err;
    Double_t v4;
    Double_t v4_err;
    VnUnit(Double_t _v2, Double_t _v2_err, Double_t _v3, Double_t _v3_err, Double_t _v4, Double_t _v4_err) :
        v2(_v2), v2_err(_v2_err), v3(_v3), v3_err(_v3_err), v4(_v4), v4_err(_v4_err) {}
};

void ProcessConfig_PtDiff(Bool_t isNch, std::vector<InputUnit> dataList, std::string outputFileName);
VnUnit* Vn3times2PC(Bool_t isNch, std::vector<InputUnit> dataList, std::vector<TFile*>& vnDeltaFiles, Double_t pTMin=0, Double_t pTMax=0);
VnUnit* GetResultsFromVnDeltaFiles(int isample, TFile* file, Bool_t isPtDiff, Int_t minRange, Int_t maxRange,  Double_t pTMin=0, Double_t pTMax=0);
double Get3times2PC(double vnDelta_LM, double vnDelta_MR, double vnDelta_LR);
double Get3times2PC_Error(double LM, double LM_err, double MR, double MR_err, double LR, double LR_err);
VnUnit* Get3times2PC(VnUnit* LMvalues, VnUnit* MRvalues, VnUnit* LRvalues);

void Process_3times2PC() {
    // 不显示窗口
    gROOT->SetBatch(kTRUE);
    std::vector<ConfigUnit> configList;

    // LM, MR, LR
    // configList.push_back(ConfigUnit(kCent, kPtDiffOn, 
    // {InputUnit("LHC25af_pass1_532068", kTPCFT0A, kTemplateFit, 0, 10), InputUnit("LHC25af_pass1_532068", kTPCFT0C, kTemplateFit, 0, 10), InputUnit("LHC25af_pass1_537547", kFT0AFT0C, kFourierFit, 0, 20)}, 
    // "LHC25af_pass1_537547"));
    // Use TemplateFit VnDelta files: LM=TPC-FT0C, MR=TPC-FT0A, LR=FT0A-FT0C
    configList.push_back(ConfigUnit(kCent, kPtDiffOn, 
    {InputUnit("LHC25af_pass2_565246", kTPCFT0C, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_565246", kTPCFT0A, kTemplateFit, 0, 20), InputUnit("LHC25af_pass2_563812", kFT0AFT0C, kTemplateFit, 0, 20)}, 
    "LHC25af_pass2_565246_kTemplateFit"));
    // TPC eta-split example: (LM) TPC_TPC eta>0.5, (MR) TPC_TPC eta<-0.5, (LR) cross eta>0.5 × eta<-0.5
    // Ensure you have produced matching VnDelta files with these suffixes upstream (TemplateFit/FourierFit).
    // configList.push_back(ConfigUnit(kCent, kPtDiffOn,
    // {InputUnit("LHC25af_pass2_565246_EtaPos", kTPCTPC, kTemplateFit, 0, 20),
    //  InputUnit("LHC25af_pass2_565246_EtaNeg", kTPCTPC, kTemplateFit, 0, 20),
    //  InputUnit("LHC25af_pass2_565246_EtaPosEtaNeg", kTPCTPC, kTemplateFit, 0, 20)},
    // "LHC25af_pass2_565246_TPCTPC_EtaSplit_kTemplateFit"));

    for (auto config : configList) {
        if (config.isPtDiff) {
            ProcessConfig_PtDiff(config.isNch, config.dataList, config.outputFileName);
        }
    }
}


void ProcessConfig_PtDiff(Bool_t isNch, std::vector<InputUnit> dataList, std::string outputFileName) {
    if (dataList.size() != 3) {
        Printf("The dataList don't have a size of 3!");
        return;
    }
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    std::vector<TFile*> vnDeltaFiles(3);
    for (int index = 0; index < 2; index++) {
        auto data = dataList[index];
        vnDeltaFiles[index] = new TFile(Form("./%s/PtDiff/VnDelta_%s_%s_%i_%i_%s.root", DihadronMethodName[data.method].c_str(), data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, DihadronCorrTypeName[data.corrType].c_str()), "READ");
        if (!vnDeltaFiles[index]) {
            Printf("can not get file!");
            return;
        }
    }
    {
        auto data = dataList[2];
        vnDeltaFiles[2] = new TFile(Form("./%s/VnDelta_%s_%s_%s.root", DihadronMethodName[data.method].c_str(), data.fileNameSuffix.c_str(), splitName.c_str(), DihadronCorrTypeName[data.corrType].c_str()), "READ");
        if (!vnDeltaFiles[2]) {
            Printf("can not get file for reference!");
            return;
        }
    }

    std::vector<VnUnit*> vnResults;
    for (Int_t iPt = 0; iPt < pTBins.size() - 1; iPt++) {
        double pTMin = pTBins[iPt];
        double pTMax = pTBins[iPt + 1];
        vnResults.push_back(Vn3times2PC(isNch, dataList, vnDeltaFiles, pTMin, pTMax));
    }

    TFile outputFile(Form("./3times2PC/Vn_%s_%s_%i_%i.root", outputFileName.c_str(), splitName.c_str(), dataList[2].minRange, dataList[2].maxRange), "RECREATE");

    // 初始化直方图
    TH1D* hV2 = new TH1D("hV2", "v_{2};p_{T};v_{2}", 
                        pTBins.size()-1, pTBins.data());
    TH1D* hV3 = new TH1D("hV3", "v_{3};p_{T};v_{3}", 
                        pTBins.size()-1, pTBins.data());
    TH1D* hV4 = new TH1D("hV4", "v_{4};p_{T};v_{4}", 
                        pTBins.size()-1, pTBins.data());

    // 填充数据
    for (size_t i = 0; i < vnResults.size(); ++i) {
        hV2->SetBinContent(i+1, vnResults[i]->v2);
        hV2->SetBinError(i+1, vnResults[i]->v2_err);
        
        hV3->SetBinContent(i+1, vnResults[i]->v3);
        hV3->SetBinError(i+1, vnResults[i]->v3_err);
        
        hV4->SetBinContent(i+1, vnResults[i]->v4);
        hV4->SetBinError(i+1, vnResults[i]->v4_err);
    }

    // 写入文件
    hV2->Write();
    hV3->Write();
    hV4->Write();

    std::cout << "Output file: " << Form("./3times2PC/Vn_%s_%s_%i_%i.root", outputFileName.c_str(), splitName.c_str(), dataList[2].minRange, dataList[2].maxRange) << std::endl;
    outputFile.Close();
}

VnUnit* Vn3times2PC(Bool_t isNch, std::vector<InputUnit> dataList, std::vector<TFile*>& vnDeltaFiles, Double_t pTMin=0, Double_t pTMax=0) {
    if (dataList.size() != 3) {
        Printf("The dataList don't have a size of 3!");
        return;
    }
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    // central values
    VnUnit* resultsLM = nullptr;
    VnUnit* resultsMR = nullptr;
    VnUnit* resultsLR = nullptr;
    resultsLM = GetResultsFromVnDeltaFiles(-1, vnDeltaFiles[0], kTRUE, dataList[0].minRange, dataList[0].maxRange, pTMin, pTMax);
    resultsMR = GetResultsFromVnDeltaFiles(-1, vnDeltaFiles[1], kTRUE, dataList[1].minRange, dataList[1].maxRange, pTMin, pTMax);
    resultsLR = GetResultsFromVnDeltaFiles(-1, vnDeltaFiles[2], kFALSE, dataList[2].minRange, dataList[2].maxRange, pTMin, pTMax);
    if (!resultsLM || !resultsMR || !resultsLR) {
        Printf("Can not get results from vnDelta files!");
        return;
    }
    VnUnit* vnResult = Get3times2PC(resultsLM, resultsMR, resultsLR);

    std::vector<std::vector<std::vector<double>>> ValueArray;
    std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    std::vector<std::vector<double>> ErrorArray;
    int Nobs=3;//v22,v32,v42
    int NofSample = maxSample*maxSample;
    int Nbin = 1;
    ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);
    for(int sample=0;sample<NofSample;sample++) {
        VnUnit* tempLM = nullptr;
        VnUnit* tempMR = nullptr;
        VnUnit* tempLR = nullptr;
        tempLM = GetResultsFromVnDeltaFiles(sample, vnDeltaFiles[0], kTRUE, dataList[0].minRange, dataList[0].maxRange, pTMin, pTMax);
        tempMR = GetResultsFromVnDeltaFiles(sample, vnDeltaFiles[1], kTRUE, dataList[1].minRange, dataList[1].maxRange, pTMin, pTMax);
        tempLR = GetResultsFromVnDeltaFiles(sample, vnDeltaFiles[2], kFALSE, dataList[2].minRange, dataList[2].maxRange, pTMin, pTMax);
        if (!tempLM || !tempMR || !tempLR) {
            Printf("Can not get subsample %d results from vnDelta files!", sample);
            return;
        }
        VnUnit* vnTemp = Get3times2PC(tempLM, tempMR, tempLR);
        ValueArray[0][sample][0] = vnTemp->v2;
        ValueErrorArray[0][sample][0] = vnTemp->v2_err;
        ValueArray[1][sample][0] = vnTemp->v3;
        ValueErrorArray[1][sample][0] = vnTemp->v3_err;
        ValueArray[2][sample][0] = vnTemp->v4;
        ValueErrorArray[2][sample][0] = vnTemp->v4_err;
        delete tempLM;
        delete tempMR;
        delete tempLR;
        delete vnTemp;
    }

    for(int iobs = 0;iobs < Nobs;iobs++){
        CalculateBootstrapError(ValueArray[iobs],ValueErrorArray[iobs],ErrorArray[iobs],1.);
    }

    vnResult->v2_err = ErrorArray[0][0];
    vnResult->v3_err = ErrorArray[1][0];
    vnResult->v4_err = ErrorArray[2][0];

    std::cout << "print result: " << dataList[2].fileNameSuffix << std::endl;
    std::cout << "v2: " << vnResult->v2 << " +/- " << vnResult->v2_err << std::endl;
    std::cout << "v3: " << vnResult->v3 << " +/- " << vnResult->v3_err << std::endl;
    std::cout << "v4: " << vnResult->v4 << " +/- " << vnResult->v4_err << std::endl;

    return vnResult;
}


VnUnit* GetResultsFromVnDeltaFiles(int isample, TFile* file, Bool_t isPtDiff, Int_t minRange, Int_t maxRange,  Double_t pTMin=0, Double_t pTMax=0) {
    TH1D* hv2 = nullptr;
    TH1D* hv3 = nullptr;
    TH1D* hv4 = nullptr;
    double v2 = -1;
    double v3 = -1;
    double v4 = -1;
    double v2_err = 10;
    double v3_err = 10;
    double v4_err = 10;
    if(isample > 0) {
        hv2 = (TH1D*)file->Get(Form("Subsamples/hV2_subsample_%d", isample));
        hv3 = (TH1D*)file->Get(Form("Subsamples/hV3_subsample_%d", isample));
        hv4 = (TH1D*)file->Get(Form("Subsamples/hV4_subsample_%d", isample));
    } else {
        hv2 = (TH1D*)file->Get("hV2");
        hv3 = (TH1D*)file->Get("hV3");
        hv4 = (TH1D*)file->Get("hV4");
    }
    if (!hv2 || !hv3 || !hv4) {
        Printf("Can not get histograms from file!");
        return;
    }
    if (isPtDiff) {
        // Get bincontent from the bin of [pTMin, pTMax] (exactly that one bin)
        int binIndex = hv2->GetXaxis()->FindBin((pTMin + pTMax) / 2.0);
        v2 = hv2->GetBinContent(binIndex);
        v2_err = hv2->GetBinError(binIndex);
        v3 = hv3->GetBinContent(binIndex);
        v3_err = hv3->GetBinError(binIndex);
        v4 = hv4->GetBinContent(binIndex);
        v4_err = hv4->GetBinError(binIndex);
    } else {
        // Get bincontent from the bin of [minRange, maxRange]
        int binIndex = hv2->GetXaxis()->FindBin((minRange + maxRange) / 2.0);
        v2 = hv2->GetBinContent(binIndex);
        v2_err = hv2->GetBinError(binIndex);
        v3 = hv3->GetBinContent(binIndex);
        v3_err = hv3->GetBinError(binIndex);
        v4 = hv4->GetBinContent(binIndex);
        v4_err = hv4->GetBinError(binIndex);
    }
    VnUnit* vnResult = new VnUnit(v2, v2_err, v3, v3_err, v4, v4_err);
    return vnResult;
}

double Get3times2PC(double vnDelta_LM, double vnDelta_MR, double vnDelta_LR) {
    double vn = (vnDelta_LM * vnDelta_MR) / vnDelta_LR;
    if (vn >= 0) {
        return std::sqrt(vn);
    } else {
        return -1;
    }
}

double Get3times2PC_Error(double LM, double LM_err, double MR, double MR_err, double LR, double LR_err) {
    if (LM < 0 || MR < 0 || LR < 0) {
        return 10.;
    }
    double part1 = 0.5*pow(LM, -0.5)*pow(MR, 0.5)*pow(LR,-0.5)*LM_err;
    part1 = pow(part1, 2);
    double part2 = 0.5*pow(LM, 0.5)*pow(MR, -0.5)*pow(LR,-0.5)*MR_err;
    part2 = pow(part2, 2);
    double part3 = 0.5*pow(LM, 0.5)*pow(MR, 0.5)*pow(LR,-1.5)*LR_err;
    part3 = pow(part3, 2);
    return sqrt(part1 + part2 + part3);
}

VnUnit* Get3times2PC(VnUnit* LMvalues, VnUnit* MRvalues, VnUnit* LRvalues) {
    double v2 = Get3times2PC(LMvalues->v2, MRvalues->v2, LRvalues->v2);
    double v2_err = Get3times2PC_Error(LMvalues->v2, LMvalues->v2_err, MRvalues->v2, MRvalues->v2_err, LRvalues->v2, LRvalues->v2_err);
    double v3 = Get3times2PC(LMvalues->v3, MRvalues->v3, LRvalues->v3);
    double v3_err = Get3times2PC_Error(LMvalues->v3, LMvalues->v3_err, MRvalues->v3, MRvalues->v3_err, LRvalues->v3, LRvalues->v3_err);
    double v4 = Get3times2PC(LMvalues->v4, MRvalues->v4, LRvalues->v4);
    double v4_err = Get3times2PC_Error(LMvalues->v4, LMvalues->v4_err, MRvalues->v4, MRvalues->v4_err, LRvalues->v4, LRvalues->v4_err);
    VnUnit* vnResult = new VnUnit(v2, v2_err, v3, v3_err, v4, v4_err);
    return vnResult;
}