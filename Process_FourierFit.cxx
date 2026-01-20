/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch)  
 * @Date: 2025-06-23 11:16:34 
 * @Last Modified by: Zhiyong Lu
 * @Last Modified time: 2025-10-24 17:03:02
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
#include "TF1.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "./include/BasicForDihadron.h"
#include "./include/Bootstrap.h"
#include "./include/plotting.h"
#include "./include/ErrorPropagation.h"

// define struct
struct InputUnit {
    std::string fileNameSuffix;
    Int_t corrType;
    Int_t minRange;
    Int_t maxRange;
    InputUnit(std::string _fileNameSuffix, Int_t _corrType, Int_t _minRange, Int_t _maxRange) :
        fileNameSuffix(_fileNameSuffix), corrType(_corrType), minRange(_minRange), maxRange(_maxRange) {}
};

struct ConfigUnit {
    Bool_t isNch;
    Bool_t isPtDiff;
    std::vector<InputUnit> dataList;
    std::string outputFileName;
    ConfigUnit(Bool_t _isNch, Bool_t _isPtDiff, std::vector<InputUnit> _dataList, std::string _outputFileName) :
        isNch(_isNch), isPtDiff(_isPtDiff), dataList(_dataList), outputFileName(_outputFileName) {}
};

struct VnUnit {
    Double_t v2;
    Double_t v2_err;
    Double_t v3;
    Double_t v3_err;
    Double_t v4;
    Double_t v4_err;
    std::vector<Double_t> subsample_v2;
    std::vector<Double_t> subsample_v2_err;
    std::vector<Double_t> subsample_v3;
    std::vector<Double_t> subsample_v3_err;
    std::vector<Double_t> subsample_v4;
    std::vector<Double_t> subsample_v4_err;
    VnUnit(Double_t _v2, Double_t _v2_err, Double_t _v3, Double_t _v3_err, Double_t _v4, Double_t _v4_err) :
        v2(_v2), v2_err(_v2_err), v3(_v3), v3_err(_v3_err), v4(_v4), v4_err(_v4_err) {}
    void ResizeSubsample(int Nsample) {
        subsample_v2.resize(Nsample);
        subsample_v2_err.resize(Nsample);
        subsample_v3.resize(Nsample);
        subsample_v3_err.resize(Nsample);
        subsample_v4.resize(Nsample);
        subsample_v4_err.resize(Nsample);
    }
    void Fillsample(int isample, Double_t _v2, Double_t _v2_err, Double_t _v3, Double_t _v3_err, Double_t _v4, Double_t _v4_err) {
        if (isample >= subsample_v2.size()) {
            Printf("You are filling out of range: isample %d, vector size: %zu", isample, subsample_v2.size());
            return;
        }
        subsample_v2[isample] = _v2;
        subsample_v2_err[isample] = _v2_err;
        subsample_v3[isample] = _v3;
        subsample_v3_err[isample] = _v3_err;
        subsample_v4[isample] = _v4;
        subsample_v4_err[isample] = _v4_err;
    }
};

// declare functions
void ProcessConfig(Bool_t isNch, std::vector<InputUnit> dataList, std::string outputFileName);
void ProcessConfig_PtDiff(Bool_t isNch, std::vector<InputUnit> dataList, std::string outputFileName);
VnUnit* FourierFit(Bool_t isNch, InputUnit data, Bool_t cn2Tovn2, Double_t pTMin=0, Double_t pTMax=0);
std::vector<Int_t> CheckAndMergeRanges(const std::vector<InputUnit>& inputUnits);
VnUnit* fitSample(Bool_t isNch, TFile* datafile, InputUnit data, int sample = -1, Double_t pTMin=0, Double_t pTMax=0);
void RooFourierFit(TH1 *hist, std::vector<Double_t>& fParamVal, std::vector<Double_t>& fParamErr);
void PlotFitting(TH1 *hm, Bool_t isNch, std::string fileSuffix, Int_t minRange, Int_t maxRange, const std::vector<Double_t>& par, const std::vector<Double_t>& parerr, Int_t corrType, Double_t pTMin=0, Double_t pTMax=0);

// global variables
std::string collisionSystemName = "peripheral PbPb";
Bool_t kOutputVnDelta = true;
//==============================================================
void Process_FourierFit() {
    // 不显示窗口
    gROOT->SetBatch(kTRUE);
    std::vector<ConfigUnit> configList;

    collisionSystemName = "Ne-Ne";
    kOutputVnDelta = true;
    //old
    configList.push_back(ConfigUnit(kCent, kPtDiffOff,
    {InputUnit("LHC25af_pass2_565246", kTPCFT0A, 0, 20)}, 
    "LHC25af_pass2_565246"));
    configList.push_back(ConfigUnit(kCent, kPtDiffOff,
    {InputUnit("LHC25af_pass2_565246", kTPCFT0C, 0, 20)}, 
    "LHC25af_pass2_565246"));
    //remap
    configList.push_back(ConfigUnit(kCent, kPtDiffOff,
    {InputUnit("LHC25af_pass2_591069", kTPCFT0A, 0, 20)}, 
    "LHC25af_pass2_591069"));
    configList.push_back(ConfigUnit(kCent, kPtDiffOff,
    {InputUnit("LHC25af_pass2_591069", kTPCFT0C, 0, 20)}, 
    "LHC25af_pass2_591069"));
    //innerring
    configList.push_back(ConfigUnit(kCent, kPtDiffOff,
    {InputUnit("LHC25af_pass2_591387", kTPCFT0A, 0, 20)}, 
    "LHC25af_pass2_591387"));
    configList.push_back(ConfigUnit(kCent, kPtDiffOff,
    {InputUnit("LHC25af_pass2_591387", kTPCFT0C, 0, 20)}, 
    "LHC25af_pass2_591387"));
    //outterring
    configList.push_back(ConfigUnit(kCent, kPtDiffOff,
    {InputUnit("LHC25af_pass2_591389", kTPCFT0A, 0, 20)}, 
    "LHC25af_pass2_591389"));
    configList.push_back(ConfigUnit(kCent, kPtDiffOff,
    {InputUnit("LHC25af_pass2_591389", kTPCFT0C, 0, 20)}, 
    "LHC25af_pass2_591389"));
    
    for (auto config : configList) {
        if (config.isPtDiff) {
            ProcessConfig_PtDiff(config.isNch, config.dataList, config.outputFileName);
        } else {
            ProcessConfig(config.isNch, config.dataList, config.outputFileName);
        }
    }

}

//==============================================================
void ProcessConfig(Bool_t isNch, std::vector<InputUnit> dataList, std::string outputFileName) {
    // 按 minRange 排序输入列表
    std::sort(dataList.begin(), dataList.end(), [](const InputUnit& a, const InputUnit& b) {
        return a.minRange < b.minRange;
    });
    
    // print datalist
    std::cout << "Data list: " << std::endl;
    for (const auto& data : dataList) {
        std::cout << "[" << data.minRange << ", " << data.maxRange << "] " << std::endl;
    }

    // 检查范围是否连续
    std::vector<Int_t> mergedRanges = CheckAndMergeRanges(dataList);
    Bool_t isContinuous = !mergedRanges.empty();

    // 执行模板拟合获取所有结果
    std::vector<VnUnit*> vnResults;
    for (const auto& data : dataList) {
        vnResults.push_back(FourierFit(isNch, data, (!kOutputVnDelta)));
    }

    // 创建输出文件
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    std::string stringDelta = "";
    if (kOutputVnDelta) stringDelta = "Delta";
    TFile outputFile(Form("./FourierFit/Vn%s_%s_%s_%s.root", stringDelta.c_str(), outputFileName.c_str(), splitName.c_str(), DihadronCorrTypeName[dataList[0].corrType].c_str()), "RECREATE");

    if (isContinuous) {
        // 创建可变bin宽度的TH1D
        std::vector<Double_t> binEdges;
        for (auto val : mergedRanges) {
            binEdges.push_back(static_cast<Double_t>(val));
        }

        // 初始化直方图
        TH1D* hV2 = new TH1D("hV2", "v_{2};Centrality;v_{2}", 
                            mergedRanges.size()-1, binEdges.data());
        TH1D* hV3 = new TH1D("hV3", "v_{3};Centrality;v_{3}", 
                            mergedRanges.size()-1, binEdges.data());
        TH1D* hV4 = new TH1D("hV4", "v_{4};Centrality;v_{4}", 
                            mergedRanges.size()-1, binEdges.data());

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
        if (kOutputVnDelta) {
            // create sub directory and write the subsample results
            TDirectory* subsampleDir = outputFile.mkdir("Subsamples");
            subsampleDir->cd();
            for (Int_t sample = 0; sample < vnResults[0]->subsample_v2.size(); sample++) {
                TH1D* hV2_sub = new TH1D(Form("hV2_subsample_%d", sample), Form("v_{2} subsample %d;Centrality;v_{2}", sample), 
                                    mergedRanges.size()-1, binEdges.data());
                TH1D* hV3_sub = new TH1D(Form("hV3_subsample_%d", sample), Form("v_{3} subsample %d;Centrality;v_{3}", sample), 
                                    mergedRanges.size()-1, binEdges.data());
                TH1D* hV4_sub = new TH1D(Form("hV4_subsample_%d", sample), Form("v_{4} subsample %d;Centrality;v_{4}", sample), 
                                    mergedRanges.size()-1, binEdges.data());
                for (size_t i = 0; i < vnResults.size(); ++i) {
                    hV2_sub->SetBinContent(i+1, vnResults[i]->subsample_v2[sample]);
                    hV2_sub->SetBinError(i+1, vnResults[i]->subsample_v2_err[sample]);
                    hV3_sub->SetBinContent(i+1, vnResults[i]->subsample_v3[sample]);
                    hV3_sub->SetBinError(i+1, vnResults[i]->subsample_v3_err[sample]);
                    hV4_sub->SetBinContent(i+1, vnResults[i]->subsample_v4[sample]);
                    hV4_sub->SetBinError(i+1, vnResults[i]->subsample_v4_err[sample]);
                }
                hV2_sub->Write();
                hV3_sub->Write();
                hV4_sub->Write();
            }
        }
    } else {
        // 创建TGraphErrors
        Int_t nPoints = dataList.size();
        TGraphErrors* gV2 = new TGraphErrors(nPoints);
        TGraphErrors* gV3 = new TGraphErrors(nPoints);
        TGraphErrors* gV4 = new TGraphErrors(nPoints);

        // 设置标题
        std::string Xtitle = "Centrality (%)";
        if (isNch) Xtitle = "Nch";
        gV2->SetNameTitle("gV2", Form("v_{2};%s;v_{2}", Xtitle.c_str()));
        gV3->SetNameTitle("gV3", Form("v_{3};%s;v_{3}", Xtitle.c_str()));
        gV4->SetNameTitle("gV4", Form("v_{4};%s;v_{4}", Xtitle.c_str()));

        // 填充数据
        for (Int_t i = 0; i < nPoints; ++i) {
            Double_t xCenter = 0.5*(dataList[i].minRange + dataList[i].maxRange);
            Double_t xError = 0.5*(dataList[i].maxRange - dataList[i].minRange);
            
            gV2->SetPoint(i, xCenter, vnResults[i]->v2);
            gV2->SetPointError(i, xError, vnResults[i]->v2_err);
            
            gV3->SetPoint(i, xCenter, vnResults[i]->v3);
            gV3->SetPointError(i, xError, vnResults[i]->v3_err);
            
            gV4->SetPoint(i, xCenter, vnResults[i]->v4);
            gV4->SetPointError(i, xError, vnResults[i]->v4_err);
        }

        // 写入文件
        gV2->Write();
        gV3->Write();
        gV4->Write();
    }

    std::cout << "Output file: " << Form("./FourierFit/Vn%s_%s_%s_%s.root", stringDelta.c_str(), outputFileName.c_str(), splitName.c_str(), DihadronCorrTypeName[dataList[0].corrType].c_str()) << std::endl;
    outputFile.Close();
}

//==============================================================
void ProcessConfig_PtDiff(Bool_t isNch, std::vector<InputUnit> dataList, std::string outputFileName) {
    // Just looping the list
    for (const auto& data : dataList) {
        
        // 执行模板拟合获取所有结果
        std::vector<VnUnit*> vnResults;
        for (Int_t iPt = 0; iPt < pTBins.size() - 1; iPt++) {
            double pTMin = pTBins[iPt];
            double pTMax = pTBins[iPt + 1];
            vnResults.push_back(FourierFit(isNch, data, kFALSE, pTMin, pTMax));
        }
        // 创建输出文件
        std::string splitName = "Mult";
        if (!isNch) splitName = "Cent";
        std::string stringDelta = "";
        if (kOutputVnDelta) stringDelta = "Delta";
        TFile outputFile(Form("./FourierFit/PtDiff/Vn%s_%s_%s_%i_%i_%s.root", stringDelta.c_str(), outputFileName.c_str(), splitName.c_str(), data.minRange, data.maxRange, DihadronCorrTypeName[dataList[0].corrType].c_str()), "RECREATE");

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

        if (kOutputVnDelta) {
            // create sub directory and write the subsample results
            TDirectory* subsampleDir = outputFile.mkdir("Subsamples");
            subsampleDir->cd();
            for (Int_t sample = 0; sample < vnResults[0]->subsample_v2.size(); sample++) {
                TH1D* hV2_sub = new TH1D(Form("hV2_subsample_%d", sample), "v_{2};p_{T};v_{2}", 
                                    pTBins.size()-1, pTBins.data());
                TH1D* hV3_sub = new TH1D(Form("hV3_subsample_%d", sample), "v_{3};p_{T};v_{3}", 
                                    pTBins.size()-1, pTBins.data());
                TH1D* hV4_sub = new TH1D(Form("hV4_subsample_%d", sample), "v_{4};p_{T};v_{4}", 
                                    pTBins.size()-1, pTBins.data());
                for (size_t i = 0; i < vnResults.size(); ++i) {
                    hV2_sub->SetBinContent(i+1, vnResults[i]->subsample_v2[sample]);
                    hV2_sub->SetBinError(i+1, vnResults[i]->subsample_v2_err[sample]);
                    hV3_sub->SetBinContent(i+1, vnResults[i]->subsample_v3[sample]);
                    hV3_sub->SetBinError(i+1, vnResults[i]->subsample_v3_err[sample]);
                    hV4_sub->SetBinContent(i+1, vnResults[i]->subsample_v4[sample]);
                    hV4_sub->SetBinError(i+1, vnResults[i]->subsample_v4_err[sample]);
                }
                hV2_sub->Write();
                hV3_sub->Write();
                hV4_sub->Write();
            }
        }

        std::cout << "Output file: " << Form("./FourierFit/PtDiff/Vn%s_%s_%s_%i_%i_%s.root", stringDelta.c_str(), outputFileName.c_str(), splitName.c_str(), data.minRange, data.maxRange, DihadronCorrTypeName[dataList[0].corrType].c_str()) << std::endl;
        outputFile.Close();
    }
    
}

void VnPtDiff(VnUnit* VnResult_PtDiff, VnUnit* VnResult_Ref) {
    if (!VnResult_PtDiff || !VnResult_Ref) {
        std::cerr << "Invalid VnUnit pointers provided." << std::endl;
        return;
    }
    // assuming V_{n\Delta}
    if (VnResult_Ref->v2 > 0) {
        double v2_2PC = VnResult_PtDiff->v2 / TMath::Sqrt(VnResult_Ref->v2);
        // synthetic error
        double v2_2PC_err = Error_Ratio_sqrtY(VnResult_PtDiff->v2, VnResult_PtDiff->v2_err, VnResult_Ref->v2, VnResult_Ref->v2_err);
        VnResult_PtDiff->v2 = v2_2PC;
        VnResult_PtDiff->v2_err = v2_2PC_err;
    } else {
        VnResult_PtDiff->v2 = -1.;
        VnResult_PtDiff->v2_err = 10.;
    }
    if (VnResult_Ref->v3 > 0) {
        double v3_2PC = VnResult_PtDiff->v3 / TMath::Sqrt(VnResult_Ref->v3);
        // synthetic error
        double v3_2PC_err = Error_Ratio_sqrtY(VnResult_PtDiff->v3, VnResult_PtDiff->v3_err, VnResult_Ref->v3, VnResult_Ref->v3_err);
        VnResult_PtDiff->v3 = v3_2PC;
        VnResult_PtDiff->v3_err = v3_2PC_err;
    } else {
        VnResult_PtDiff->v3 = -1.;
        VnResult_PtDiff->v3_err = 10.;
    }
    if (VnResult_Ref->v4 > 0) {
        double v4_2PC = VnResult_PtDiff->v4 / TMath::Sqrt(VnResult_Ref->v4);
        // synthetic error
        double v4_2PC_err = Error_Ratio_sqrtY(VnResult_PtDiff->v4, VnResult_PtDiff->v4_err, VnResult_Ref->v4, VnResult_Ref->v4_err);
        VnResult_PtDiff->v4 = v4_2PC;
        VnResult_PtDiff->v4_err = v4_2PC_err;
    } else {
        VnResult_PtDiff->v4 = -1.;
        VnResult_PtDiff->v4_err = 10.;
    }
    return;
}

//==============================================================
VnUnit* FourierFit(Bool_t isNch, InputUnit data, Bool_t cn2Tovn2, Double_t pTMin=0, Double_t pTMax=0) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    TFile* datafile = new TFile(Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d_%s.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, DihadronCorrTypeName[data.corrType].c_str()), "READ");
    if (!datafile || !datafile->IsOpen()) {
        std::cerr << "Cannot open input file: " << Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d_%s.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, DihadronCorrTypeName[data.corrType].c_str()) << std::endl;
        exit(1);
    }
    TFile* datafile_PtDiff = nullptr;
    if (pTMin > 0 && pTMax > 0) {
        datafile_PtDiff = new TFile(Form("./ProcessOutput/PtDiff/BootstrapSample_%s_%s_%i_%i_Pt_%0.1f_%0.1f_%s.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, pTMin, pTMax, DihadronCorrTypeName[data.corrType].c_str()), "READ");
        if (!datafile_PtDiff || !datafile_PtDiff->IsOpen()) {
            std::cerr << "Cannot open input file: " << Form("./ProcessOutput/PtDiff/BootstrapSample_%s_%s_%i_%i_Pt_%0.1f_%0.1f_%s.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, pTMin, pTMax, DihadronCorrTypeName[data.corrType].c_str()) << std::endl;
            exit(1);
        }
    }
    
    VnUnit* vnResult = fitSample(isNch, datafile, data, -1);
    if (!vnResult) {
        std::cerr << "Cannot fit sample: " << data.fileNameSuffix << std::endl;
        exit(1);
    }
    VnUnit* vnResult_PtDiff = nullptr;
    if (pTMin > 0 && pTMax > 0) {
        vnResult_PtDiff = fitSample(isNch, datafile_PtDiff, data, -1, pTMin, pTMax);
        if (!vnResult_PtDiff) {
            std::cerr << "Cannot fit pT-diff sample: " << data.fileNameSuffix << std::endl;
            exit(1);
        }
        if (!kOutputVnDelta) {
            VnPtDiff(vnResult_PtDiff, vnResult);
        }
    }

    std::vector<std::vector<std::vector<double>>> ValueArray;
    std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    std::vector<std::vector<double>> ErrorArray;
    int Nobs=3;//v22,v32,v42
    int NofSample = maxSample*maxSample;
    int Nbin = 1;
    ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);
    vnResult->ResizeSubsample(NofSample);
    if (pTMin > 0 && pTMax > 0) vnResult_PtDiff->ResizeSubsample(NofSample);

    for(int sample=0;sample<NofSample;sample++) {
        VnUnit* vnTemp = fitSample(isNch, datafile, data, sample);
        if (!vnTemp) {
            std::cerr << "Cannot fit sample: " << data.fileNameSuffix << " sample: " << sample << std::endl;
            exit(1);
        }
        ValueArray[0][sample][0] = vnTemp->v2;
        ValueErrorArray[0][sample][0] = vnTemp->v2_err;
        ValueArray[1][sample][0] = vnTemp->v3;
        ValueErrorArray[1][sample][0] = vnTemp->v3_err;
        ValueArray[2][sample][0] = vnTemp->v4;
        ValueErrorArray[2][sample][0] = vnTemp->v4_err;
        vnResult->Fillsample(sample, vnTemp->v2, vnTemp->v2_err, vnTemp->v3, vnTemp->v3_err, vnTemp->v4, vnTemp->v4_err);
        if (pTMin > 0 && pTMax > 0) {
            VnUnit* vnTemp_PtDiff = fitSample(isNch, datafile_PtDiff, data, sample, pTMin, pTMax);
            if (!vnTemp_PtDiff) {
                std::cerr << "Cannot fit pT-diff sample: " << data.fileNameSuffix << std::endl;
                exit(1);
            }
            if (!kOutputVnDelta) {
                VnPtDiff(vnTemp_PtDiff, vnTemp);
            }
            ValueArray[0][sample][0] = vnTemp_PtDiff->v2;
            ValueErrorArray[0][sample][0] = vnTemp_PtDiff->v2_err;
            ValueArray[1][sample][0] = vnTemp_PtDiff->v3;
            ValueErrorArray[1][sample][0] = vnTemp_PtDiff->v3_err;
            ValueArray[2][sample][0] = vnTemp_PtDiff->v4;
            ValueErrorArray[2][sample][0] = vnTemp_PtDiff->v4_err;
            vnResult_PtDiff->Fillsample(sample, vnTemp_PtDiff->v2, vnTemp_PtDiff->v2_err, vnTemp_PtDiff->v3, vnTemp_PtDiff->v3_err, vnTemp_PtDiff->v4, vnTemp_PtDiff->v4_err);
            delete vnTemp_PtDiff;
        }
        delete vnTemp;
    }
    // for(int sample=0;sample<NofSample;sample++) {
    //     std::cout << "sample: " << sample << " v2^2: " << ValueArray[0][sample][0]  << std::endl;
    //     std::cout << "sample: " << sample << " v3^2: " << ValueArray[1][sample][0]  << std::endl;
    //     std::cout << "sample: " << sample << " v4^2: " << ValueArray[2][sample][0] << std::endl << std::endl;
    // }
    for(int iobs = 0;iobs < Nobs;iobs++){
        CalculateBootstrapError(ValueArray[iobs],ValueErrorArray[iobs],ErrorArray[iobs],1.);
    }

    vnResult->v2_err = ErrorArray[0][0];
    vnResult->v3_err = ErrorArray[1][0];
    vnResult->v4_err = ErrorArray[2][0];
    if (pTMin > 0 && pTMax > 0) {
        vnResult_PtDiff->v2_err = ErrorArray[0][0];
        vnResult_PtDiff->v3_err = ErrorArray[1][0];
        vnResult_PtDiff->v4_err = ErrorArray[2][0];
    }

    if (cn2Tovn2) {
        if (vnResult->v2 > 0.) {
            vnResult->v2_err = vnResult->v2_err / (2 * sqrt(vnResult->v2));
            vnResult->v2 = sqrt(vnResult->v2);
        }
        else {
            vnResult->v2 = -1;
            vnResult->v2_err = 10.;
        }
        
        if (vnResult->v3 > 0.) {
            vnResult->v3_err = vnResult->v3_err / (2 * sqrt(vnResult->v3));
            vnResult->v3 = sqrt(vnResult->v3);
        }
        else {
            vnResult->v3 = -1;
            vnResult->v3_err = 10.;
        }
        
        if (vnResult->v4 > 0.) {
            vnResult->v4_err = vnResult->v4_err / (2 * sqrt(vnResult->v4));
            vnResult->v4 = sqrt(vnResult->v4);
        }
        else {
            vnResult->v4 = -1;
            vnResult->v4_err = 10.;
        }
    }

    // print result
    if (pTMin > 0 && pTMax > 0) {
        std::cout << "print result: " << data.fileNameSuffix << " pT-diff" << std::endl;
        std::cout << "v2: " << vnResult_PtDiff->v2 << " +/- " << vnResult_PtDiff->v2_err << std::endl;
        std::cout << "v3: " << vnResult_PtDiff->v3 << " +/- " << vnResult_PtDiff->v3_err << std::endl;
        std::cout << "v4: " << vnResult_PtDiff->v4 << " +/- " << vnResult_PtDiff->v4_err << std::endl;
        return vnResult_PtDiff;
    }
    std::cout << "print result: " << data.fileNameSuffix << std::endl;
    std::cout << "v2: " << vnResult->v2 << " +/- " << vnResult->v2_err << std::endl;
    std::cout << "v3: " << vnResult->v3 << " +/- " << vnResult->v3_err << std::endl;
    std::cout << "v4: " << vnResult->v4 << " +/- " << vnResult->v4_err << std::endl;

    return vnResult;
}


//==============================================================
std::vector<Int_t> CheckAndMergeRanges(const std::vector<InputUnit>& inputUnits) {
    std::vector<InputUnit> sortedUnits = inputUnits;
    // 按 minRange 升序排序
    std::sort(sortedUnits.begin(), sortedUnits.end(),
              [](const InputUnit& a, const InputUnit& b) {
                  return a.minRange < b.minRange;
              });

    // 检查连续性
    bool isContinuous = true;
    for (size_t i = 0; i < sortedUnits.size() - 1; ++i) {
        if (sortedUnits[i].maxRange != sortedUnits[i + 1].minRange) {
            isContinuous = false;
            break;
        }
    }

    // 如果连续则生成合并后的范围
    std::vector<Int_t> mergedRanges;
    if (isContinuous && !sortedUnits.empty()) {
        mergedRanges.reserve(sortedUnits.size() + 1);
        mergedRanges.push_back(sortedUnits[0].minRange);
        for (const auto& unit : sortedUnits) {
            mergedRanges.push_back(unit.maxRange);
        }
    }
    return mergedRanges;
}

//==============================================================
VnUnit* fitSample(Bool_t isNch, TFile* datafile, InputUnit data, int sample = -1, Double_t pTMin=0, Double_t pTMax=0) {
    std::vector<Double_t> fParamVal;
    std::vector<Double_t> fParamErr;
    TH1D* lm=0;
    TH1D* hm=0;
    TString suffix = (sample == -1) ? "" : Form("_%d", sample);
    hm = (TH1D*)datafile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d%s", data.minRange, data.maxRange, suffix.Data()));
    if (!hm) {
        std::cerr << "Cannot find histogram: " << Form("bsSample_hPhiSameOverMixed_%d_%d%s", data.minRange, data.maxRange, suffix.Data()) << std::endl;
        return 0;
    }
    RooFourierFit(hm, fParamVal, fParamErr);
    if (sample == -1) {
        PlotFitting(hm, isNch, data.fileNameSuffix, data.minRange, data.maxRange, fParamVal, fParamErr, data.corrType, pTMin, pTMax);
    }
    VnUnit* vnResult = new VnUnit(fParamVal[0], fParamErr[0], fParamVal[1], fParamErr[1], fParamVal[2], fParamErr[2]);
    return vnResult;
}

//==============================================================
void RooFourierFit(TH1 *hist, std::vector<Double_t>& fParamVal, std::vector<Double_t>& fParamErr) {
    // 检查输入直方图是否有效
    if (!hist) {
        Error("FourierFit", "Invalid histogram pointer!");
        return;
    }
    fParamVal.clear();
    fParamErr.clear();
    fParamVal.resize(5);
    fParamErr.resize(5);

    // 定义傅里叶级数函数 - 5参数形式
    auto fourierFunc = [](double *x, double *p) {
        return p[0] + 2*p[1]*TMath::Cos(x[0]) 
                    + 2*p[2]*TMath::Cos(2*x[0]) 
                    + 2*p[3]*TMath::Cos(3*x[0])
                    + 2*p[4]*TMath::Cos(4*x[0]);
    };

    // 创建TF1对象
    TF1 *fitFunc = new TF1("fourierFit", fourierFunc, 
                          hist->GetXaxis()->GetXmin(),  // 使用直方图X范围
                          hist->GetXaxis()->GetXmax(), 
                          5);  // 4个参数

    std::vector<double> Coeff;
    HistFFT(hist, Coeff);
    std::cout << "Coeff: " << Coeff[0] << " " << Coeff[1] << " " << Coeff[2] << " " << Coeff[3] << " " << Coeff[4] << " " << Coeff[5] << std::endl;
    // 设置初始参数（基于直方图特性）
    // Cent
    // fitFunc->SetParameter(0, hist->GetMean());  // 常数项初始化为平均值
    // for (int i = 1; i < 5; ++i) {
    //     fitFunc->SetParameter(i, 0.1 * hist->GetRMS());  // 谐波项初始值设为RMS的10%
    // }
    // Mult
    fitFunc->SetParameter(0, Coeff[0]);  // 常数项初始化为平均值
    for (int i = 1; i < 5; ++i) {
        fitFunc->SetParameter(i, Coeff[i]);  // 谐波项初始值设为RMS的10%
    }

    // 设置参数名称（可选）
    fitFunc->SetParName(0, "a0");
    fitFunc->SetParName(1, "a1");
    fitFunc->SetParName(2, "a2");
    fitFunc->SetParName(3, "a3");
    fitFunc->SetParName(4, "a4");

    // 执行拟合（使用最小二乘法并抑制输出）
    // hist->Fit(fitFunc, "QN");
    TFitResultPtr fitResult = hist->Fit(fitFunc, "QN0");
    // 检查拟合状态
    if (static_cast<int>(fitResult)) {
        Error("FourierFit", "Fit failed with status %d", static_cast<int>(fitResult));
        for (int i = 0; i < 5; ++i) {
            fParamVal[i] = -1;
            fParamErr[i] = 10;
        }
        return;
    }
    std::cout << "Fit results: " << std::endl;

    // 计算拟合质量指标
    double chi2 = fitFunc->GetChisquare();
    double ndf = fitFunc->GetNDF();
    double chi2_ndf = (ndf > 0) ? chi2 / ndf : 0;
    std::cout << "Chi2/NDF = " << chi2_ndf << std::endl;
    double pvalue = TMath::Prob(chi2, static_cast<int>(ndf));
    std::cout << "P-value = " << pvalue << std::endl;
    if (chi2_ndf > 2000) {
        std::cout << "WARNING: Chi2/NDF > 20, fit may be inaccurate!" << std::endl;
        for (int i = 0; i < 5; ++i) {
            fParamVal[i] = -1;
            fParamErr[i] = 10;
        }
        return;
    }


    // 获取拟合结果
    double a0 = fitFunc->GetParameter(0);
    double a0e = fitFunc->GetParError(0);
    double a1 = fitFunc->GetParameter(1);
    double a1e = fitFunc->GetParError(1);
    double a2 = fitFunc->GetParameter(2);
    double a2e = fitFunc->GetParError(2);
    double a3 = fitFunc->GetParameter(3);
    double a3e = fitFunc->GetParError(3);
    double a4 = fitFunc->GetParameter(4);
    double a4e = fitFunc->GetParError(4);
    // std::cout << "a0 = " << a0 << " +/- " << a0e << std::endl;
    // std::cout << "a1 = " << a1 << " +/- " << a1e << std::endl;
    // std::cout << "a2 = " << a2 << " +/- " << a2e << std::endl;
    // std::cout << "a3 = " << a3 << " +/- " << a3e << std::endl;
    // std::cout << "a4 = " << a4 << " +/- " << a4e << std::endl;

    double v12square = a1 / a0;
    double v12square_err = Error_Ratio(a1, a1e, a0, a0e, 0);
    double v22square = a2 / a0;
    double v22square_err = Error_Ratio(a2, a2e, a0, a0e, 0);
    double v32square = a3 / a0;
    double v32square_err = Error_Ratio(a3, a3e, a0, a0e, 0);
    double v42square = a4 / a0;
    double v42square_err = Error_Ratio(a4, a4e, a0, a0e, 0);

    std::cout << "v12^2 = " << v12square << " +/- " << v12square_err << std::endl;
    std::cout << "v22^2 = " << v22square << " +/- " << v22square_err << std::endl;
    std::cout << "v32^2 = " << v32square << " +/- " << v32square_err << std::endl;
    std::cout << "v42^2 = " << v42square << " +/- " << v42square_err << std::endl;

    
    fParamVal[0] = v22square;
    fParamErr[0] = v22square_err;
    fParamVal[1] = v32square;
    fParamErr[1] = v32square_err;
    fParamVal[2] = v42square;
    fParamErr[2] = v42square_err;
    fParamVal[3] = v12square;
    fParamErr[3] = v12square_err;
    fParamVal[4] = a0;
    fParamErr[4] = a0e;
}

void DrawText(double xmin, double ymin, double textSize, TString text)
{

	TLatex *textPreliminary = new TLatex(xmin, ymin, Form("%s", text.Data()));
	textPreliminary->SetNDC();
	textPreliminary->SetTextFont(43);
	textPreliminary->SetTextSize(textSize);
	textPreliminary->Draw();

}

void PlotFitting(TH1 *hm, Bool_t isNch, std::string fileSuffix, Int_t minRange, Int_t maxRange, const std::vector<Double_t>& par, const std::vector<Double_t>& parerr, Int_t corrType, Double_t pTMin=0, Double_t pTMax=0) {
    gStyle->SetOptStat(0); 
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    double v21 = par[0];
    double v31 = par[1];
    double v41 = par[2];
    double v12 = par[3];
    double a0 = par[4];

    double v21e = parerr[0];
    double v31e = parerr[1];
    double v41e = parerr[2];


    TCanvas* canvas = new TCanvas(Form("Fit"), "Fit", 800, 600);
    canvas->Range(0,0,1,1);
    
    // 创建上下面板
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
    pad1->SetBorderMode(0);
    pad1->SetBorderSize(2);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.05);
    pad1->SetTopMargin(0.05);
    pad1->SetBottomMargin(0.13);
    pad1->SetTicks(1,1);
    pad1->Draw();
    
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
    pad2->SetBorderMode(0);
    pad2->SetBorderSize(2);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.05);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.25);
    pad2->SetTicks(1,1);
    pad2->Draw();

    // 绘制主图
    pad1->cd();
    
    // 创建背景直方图
    TH1D* hbkg1 = new TH1D("dPhi", "dPhi", 1, -TMath::Pi()/2.0, 3*TMath::Pi()/2.0);
    hbkg1->SetStats(0);
    hbkg1->GetXaxis()->SetTitle("#Delta#phi [rad]");
    hbkg1->GetYaxis()->SetTitle("Arbitrary scaled");
    hbkg1->GetXaxis()->SetTitleSize(0.045);
    hbkg1->GetXaxis()->SetLabelSize(0.04);
    hbkg1->GetYaxis()->SetTitleSize(0.045);
    hbkg1->GetYaxis()->SetLabelSize(0.04);
    hbkg1->GetXaxis()->SetTitleOffset(1.2);
    hbkg1->GetYaxis()->SetTitleOffset(1.2);
    
    // 设置Y轴范围
    double ymax = hm->GetMaximum();
    double ymin = hm->GetMinimum();
    double ydelta = ymax - ymin;
    hbkg1->GetYaxis()->SetRangeUser(ymin-ydelta*0.2, ymax+ydelta*0.7);
    hbkg1->Draw();

    // 绘制原始数据
    hm->SetMarkerStyle(20);
    hm->SetMarkerColor(kBlack);
    hm->SetLineColor(kBlack);
    hm->SetMarkerSize(1.0);
    hm->Draw("same p");

    double a1 = v12 * a0;
    double a2 = v21 * a0;
    double a3 = v31 * a0;
    double a4 = v41 * a0;

    // =============== 新增：定义总拟合函数 ===============
    TF1* fit_p1m = new TF1("fit_p1m","[0] + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    fit_p1m->SetParameters(a0, a1, a2, a3, a4);
    
    // =============== 新增：计算chi2/ndf ===============
    double chi2 = 0.0;
    int nBins = hm->GetNbinsX();
    int nParams = 5; // 参数个数: a0, a1, a2, a3, a4
    int ndf = nBins - nParams;
    
    for (int i = 1; i <= nBins; i++) {
        double data = hm->GetBinContent(i);
        double error = hm->GetBinError(i);
        double x = hm->GetBinCenter(i);
        double fit = fit_p1m->Eval(x);
        
        if (error > 0) { // 忽略误差为0的bin
            double residual = data - fit;
            chi2 += (residual * residual) / (error * error);
        }
    }
    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

    // 创建拟合曲线
    const int pointBin = (int)hm->GetNbinsX();
    Double_t CopyPointX[pointBin];
    Double_t CopyPointY[pointBin];
    for (int i=0; i<pointBin; ++i){
      CopyPointX[i] = hm->GetBinCenter(i+1);
      double x = hm->GetBinCenter(i+1);
      CopyPointY[i] = a0 + 2*a1*cos(x) + 2*a2*cos(2*x) + 2*a3*cos(3*x) + 2*a4*cos(4*x);
    };
    TGraph* gCopy = new TGraph(pointBin,CopyPointX,CopyPointY);

    Double_t PeriPointY[pointBin];
    for (int i=0; i<pointBin; ++i){
      PeriPointY[i] = a0;
    };
    TGraph* gPeri = new TGraph(pointBin,CopyPointX,PeriPointY);

    gCopy->SetLineColor(colors[0]);
    gCopy->SetLineWidth(2);
    gCopy->Draw("same");

    gPeri->SetLineColor(colors[1]);
    gPeri->SetLineWidth(2);
    gPeri->SetLineStyle(5);
    gPeri->Draw("same");


    TF1* fit_p1 = new TF1("fit_p1","[0] + 2*[1]*cos(x)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p2 = new TF1("fit_p2","[0] + 2*[1]*cos(2*x)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p3 = new TF1("fit_p3","[0] + 2*[1]*cos(3*x)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p4 = new TF1("fit_p4","[0] + 2*[1]*cos(4*x)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    fit_p1->SetParameters(a0, a1);
    fit_p2->SetParameters(a0, a2);
    fit_p3->SetParameters(a0, a3);
    fit_p4->SetParameters(a0, a4);
    fit_p1->SetLineColor(colors[1]);
    fit_p1->SetLineWidth(2);
    fit_p1->SetLineStyle(1);
    fit_p2->SetLineColor(colors[2]);
    fit_p2->SetLineWidth(2);
    fit_p2->SetLineStyle(2);
    fit_p3->SetLineColor(colors[3]);
    fit_p3->SetLineWidth(2);
    fit_p3->SetLineStyle(3);
    fit_p4->SetLineColor(colors[4]);
    fit_p4->SetLineWidth(2);
    fit_p4->SetLineStyle(4);
    fit_p2->Draw("same");
    fit_p3->Draw("same");
    fit_p4->Draw("same");

    // 添加图例
    TLegend* leg = new TLegend(0.5, 0.65, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->AddEntry(hm, "Data", "lep");
    leg->AddEntry(gCopy, "a_{0} + #Sigma_{n=1}^{4}2a_{n}cos(n#Delta#phi)", "l");
    leg->AddEntry(gPeri, "Baseline a0", "l");
    leg->AddEntry(fit_p1, Form("a_{0} + 2a_{1}cos(#Delta#phi), v1^{2}#times10^{3} = %0.2f", v12*1e3), "l");
    leg->AddEntry(fit_p2, Form("a_{0} + 2a_{2}cos(2#Delta#phi), v2^{2}#times10^{3} = %0.2f", v21*1e3), "l");
    leg->AddEntry(fit_p3, Form("a_{0} + 2a_{3}cos(3#Delta#phi), v3^{2}#times10^{3} = %0.2f", v31*1e3), "l");
    leg->AddEntry(fit_p4, Form("a_{0} + 2a_{4}cos(4#Delta#phi), v4^{2}#times10^{3} = %0.2f", v41*1e3), "l");
    leg->Draw();

    // =============== 新增：添加chi2/ndf标签 ===============
    TLatex* chi2Label = new TLatex();
    chi2Label->SetNDC();
    chi2Label->SetTextFont(43);
    chi2Label->SetTextSize(20);
    chi2Label->DrawLatex(0.50, 0.60, Form("#chi^{2}/ndf = %.1f/%d = %.2f", chi2, ndf, chi2ndf));

    // 添加文本标签
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(43);
    tex->SetTextSize(20);
    DrawText(0.2, 0.85, 20, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
    DrawText(0.2, 0.80, 20, Form("%.1f < |#Delta#eta| < %.1f", DihadrondEtaRange[corrType][0], DihadrondEtaRange[corrType][1]));
    if (isNch) {
        DrawText(0.2, 0.75, 20, Form("%d < N_{ch} < %d", minRange, maxRange));
    }
    else {
        DrawText(0.2, 0.75, 20, Form("%d < Cent < %d", minRange, maxRange));
    }

    if (pTMin > 0 && pTMax > 0)
        DrawText(0.2, 0.70, 20, Form("%.1f < p_{T} < %.1f", pTMin, pTMax));
    DrawText(0.2, 0.65, 20, Form("V_{2#Delta} = %.5f #pm %.5f", v21, v21e));
    DrawText(0.2, 0.60, 20, Form("V_{3#Delta} = %.5f #pm %.5f", v31, v31e));

    // 绘制底部残差图
    pad2->cd();

    // TF1* fit_p1m = new TF1("fit_p1m","[0] + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    // fit_p1m->SetParameters(a0, a1, a2, a3, a4);
    
    // 创建残差直方图
    TH1D* hsubtract = (TH1D*)hm->Clone(Form("subtract"));
    // hsubtract->Add(gPeri, -1);

    TH1D* hResidual = (TH1D*)hsubtract->Clone(Form("pull"));
    double ymax_pull = 1, ymin_pull = 1;
    for (int i = 1; i <= hResidual->GetXaxis()->GetNbins(); i++) {
      double dat = hResidual->GetBinContent(i);
      double err = hResidual->GetBinError(i);
      // double lin = gCopy->GetPointY(i-1);
      double lin = fit_p1m->Eval(hResidual->GetBinCenter(i));
      hResidual->SetBinContent(i, dat/lin);
      hResidual->SetBinError(i, err/lin);
      ymin_pull = ymin_pull>(dat/lin)?(dat/lin):ymin_pull;
      ymax_pull = ymax_pull<(dat/lin)?(dat/lin):ymax_pull;
    }
    ymin_pull = 1 - (1-ymin_pull) * 2;
    ymax_pull = 1 + (ymax_pull-1) * 2;
    hResidual->GetYaxis()->SetRangeUser(ymin_pull,ymax_pull);

    // 配置残差图
    hResidual->SetMarkerStyle(20);
    hResidual->SetMarkerColor(kBlack);
    hResidual->SetLineColor(kBlack);
    hResidual->GetYaxis()->SetTitle("Data/Fit");
    hResidual->GetYaxis()->SetTitleSize(0.12);
    hResidual->GetYaxis()->SetLabelSize(0.10);
    hResidual->GetYaxis()->SetTitleOffset(0.35);
    hResidual->GetXaxis()->SetTitleSize(0.12);
    hResidual->GetXaxis()->SetTitle("#Delta#phi [rad]");
    hResidual->GetXaxis()->SetLabelSize(0.10);
    hResidual->GetXaxis()->SetTitleOffset(0.9);
    hResidual->Draw("ep");

    // 添加参考线
    TLine* line = new TLine(-TMath::Pi()/2.0, 1.0, 1.5*TMath::Pi(), 1.0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(kDashed);
    line->Draw();

    // 保存结果
    if (pTMin > 0 && pTMax > 0) {
        canvas->SaveAs(Form("./FourierFit/PtDiff/PDFs/FourierFit_%s_%s_%d_%d_Pt_%0.1f_%0.1f_%s.pdf", fileSuffix.c_str(), splitName.c_str(), minRange, maxRange, pTMin, pTMax, DihadronCorrTypeName[corrType].c_str()));
    } else {
        canvas->SaveAs(Form("./FourierFit/PDFs/FourierFit_%s_%s_%d_%d_%s.pdf", fileSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    }
    
    // 清理内存
    delete hbkg1;
    delete gCopy;
    delete gPeri;
    delete fit_p2;
    delete fit_p3;
    delete fit_p4;
    delete leg;
    delete tex;
    delete hResidual;
    delete line;
  
}
