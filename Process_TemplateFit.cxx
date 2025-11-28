/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch)  
 * @Date: 2025-05-18 13:51:01 
 * @Last Modified by: Zhiyong Lu
 * @Last Modified time: 2025-10-24 15:04:55
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
#include "./include/TemplateFitter.cxx"
#include "./include/TemplateFunction.C"
#include "./include/ErrorPropagation.h"
#include "./include/plotting.h"

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
    InputUnit templ;
    std::vector<InputUnit> dataList;
    std::string outputFileName;
    ConfigUnit(Bool_t _isNch, Bool_t _isPtDiff, InputUnit _template, std::vector<InputUnit> _dataList, std::string _outputFileName) :
        isNch(_isNch), isPtDiff(_isPtDiff), templ(_template), dataList(_dataList), outputFileName(_outputFileName) {}
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
void ProcessConfig(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName);
void ProcessConfig_PtDiff(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName);
VnUnit* TemplateFit(Bool_t isNch, InputUnit templ, InputUnit data, Bool_t cn2Tovn2, Double_t pTMin=0, Double_t pTMax=0);
std::vector<Int_t> CheckAndMergeRanges(const std::vector<InputUnit>& inputUnits);
VnUnit* fitSample(Bool_t isNch, TFile* templatefile, InputUnit templ, TFile* datafile, InputUnit data, int sample = -1, Double_t pTMin=0, Double_t pTMax=0);
void RooTempFitter(TH1 *lm, TH1 *hm, std::vector<Double_t>& fParamVal, std::vector<Double_t>& fParamErr, Bool_t kRefit);
void PlotFitting(TH1 *lm, TH1 *hm, Bool_t isNch, std::string fileSuffix, Int_t minRange, Int_t maxRange, const std::vector<Double_t>& par, const std::vector<Double_t>& parerr, Int_t corrType, Double_t pTMin=0, Double_t pTMax=0);

// global variables
std::string collisionSystemName = "peripheral PbPb";
Bool_t kOutputVnDelta = true;
//==============================================================
void Process_TemplateFit() {
    // 不显示窗口
    // gROOT->SetBatch(kTRUE);
    std::vector<ConfigUnit> configList;

    collisionSystemName = "Ne-Ne";
    kOutputVnDelta = true;
    // Ne-Ne 
    configList.push_back(ConfigUnit(kCent, kPtDiffOff, InputUnit("LHC25af_pass2_554625", kFT0AFT0C, 80, 100), 
    {InputUnit("LHC25af_pass2_556742", kFT0AFT0C, 0, 5)}, 
    "LHC25af_pass2_554625"));
    configList.push_back(ConfigUnit(kCent, kPtDiffOff, InputUnit("LHC25af_pass2_554625", kFT0AFT0C, 80, 100), 
    {InputUnit("LHC25af_pass2_556742", kFT0AFT0C, 0, 5)}, 
    "LHC25af_pass2_554625"));
    configList.push_back(ConfigUnit(kCent, kPtDiffOff, InputUnit("LHC25af_pass2_554625", kFT0AFT0C, 80, 100), 
    {InputUnit("LHC25af_pass2_563812", kFT0AFT0C, 0, 20)}, 
    "LHC25af_pass2_554625"));
    configList.push_back(ConfigUnit(kCent, kPtDiffOff, InputUnit("LHC25af_pass2_554625", kFT0AFT0C, 80, 100), 
    {InputUnit("LHC25af_pass2_563812", kFT0AFT0C, 0, 20)}, 
    "LHC25af_pass2_554625"));
    // O-O
    configList.push_back(ConfigUnit(kCent, kPtDiffOff, InputUnit("LHC25ae_pass2_561907", kFT0AFT0C, 80, 100), 
    {InputUnit("LHC25ae_pass2_561907", kFT0AFT0C, 0, 20)}, 
    "LHC25ae_pass2_561907"));
    configList.push_back(ConfigUnit(kCent, kPtDiffOff, InputUnit("LHC25ae_pass2_561907", kFT0AFT0C, 80, 100), 
    {InputUnit("LHC25ae_pass2_561907", kFT0AFT0C, 0, 20)}, 
    "LHC25ae_pass2_561907"));
    // configList.push_back(ConfigUnit(kCent, kPtDiffOn, InputUnit("LHC25af_pass1_537548", kFT0AFT0C, 80, 100), 
    // {InputUnit("LHC25af_pass1_537547", kFT0AFT0C, 0, 20)}, 
    // "LHC25af_pass1_537548"));


    for (auto config : configList) {
        if (config.isPtDiff) {
            ProcessConfig_PtDiff(config.isNch, config.templ, config.dataList, config.outputFileName);
        } else {
            ProcessConfig(config.isNch, config.templ, config.dataList, config.outputFileName);
        }
    }

}

//==============================================================
void ProcessConfig(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName) {
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
        vnResults.push_back(TemplateFit(isNch, templ, data, (!kOutputVnDelta)));
    }

    // 创建输出文件
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    std::string stringDelta = "";
    if (kOutputVnDelta) stringDelta = "Delta";
    TFile outputFile(Form("./TemplateFit/Vn%s_%s_%s_%s.root", stringDelta.c_str(), outputFileName.c_str(), splitName.c_str(), DihadronCorrTypeName[templ.corrType].c_str()), "RECREATE");

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

    std::cout << "Output file: " << Form("./TemplateFit/Vn%s_%s_%s_%s.root", stringDelta.c_str(), outputFileName.c_str(), splitName.c_str(), DihadronCorrTypeName[templ.corrType].c_str()) << std::endl;
    outputFile.Close();
}

//==============================================================
void ProcessConfig_PtDiff(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName) {
    // Just looping the list
    for (const auto& data : dataList) {
        // 执行模板拟合获取所有结果
        std::vector<VnUnit*> vnResults;
        for (Int_t iPt = 0; iPt < pTBins.size() - 1; iPt++) {
            double pTMin = pTBins[iPt];
            double pTMax = pTBins[iPt + 1];
            vnResults.push_back(TemplateFit(isNch, templ, data, kFALSE, pTMin, pTMax));
        }

        // 创建输出文件
        std::string splitName = "Mult";
        if (!isNch) splitName = "Cent";
        std::string stringDelta = "";
        if (kOutputVnDelta) stringDelta = "Delta";
        TFile outputFile(Form("./TemplateFit/PtDiff/Vn%s_%s_%s_%i_%i_%s.root", stringDelta.c_str(), outputFileName.c_str(), splitName.c_str(), data.minRange, data.maxRange, DihadronCorrTypeName[templ.corrType].c_str()), "RECREATE");

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

        std::cout << "Output file: " << Form("./TemplateFit/PtDiff/Vn%s_%s_%s_%i_%i_%s.root", stringDelta.c_str(), outputFileName.c_str(), splitName.c_str(), data.minRange, data.maxRange, DihadronCorrTypeName[templ.corrType].c_str()) << std::endl;
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
VnUnit* TemplateFit(Bool_t isNch, InputUnit templ, InputUnit data, Bool_t cn2Tovn2, Double_t pTMin=0, Double_t pTMax=0) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    // Set collisionSystemName for plotting based on dataset prefix (minimal, local change).
    // OLD: collisionSystemName was set globally at Process start and not updated per-dataset
    if (templ.fileNameSuffix.find("LHC25ae") == 0 || data.fileNameSuffix.find("LHC25ae") == 0) {
        collisionSystemName = "O-O";
    } else if (templ.fileNameSuffix.find("LHC25af") == 0 || data.fileNameSuffix.find("LHC25af") == 0) {
        collisionSystemName = "Ne-Ne";
    }

    TFile* templatefile = new TFile(Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d_%s.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange, DihadronCorrTypeName[templ.corrType].c_str()), "READ");
    if (!templatefile || !templatefile->IsOpen()) {
        std::cerr << "Cannot open template file: " << Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d_%s.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange, DihadronCorrTypeName[templ.corrType].c_str()) << std::endl;
        exit(1);
    }
    TFile* templatefile_PtDiff = nullptr;
    if (pTMin > 0 && pTMax > 0) {
        templatefile_PtDiff = new TFile(Form("./ProcessOutput/PtDiff/BootstrapSample_%s_%s_%d_%d_Pt_%0.1f_%0.1f_%s.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange, pTMin, pTMax, DihadronCorrTypeName[templ.corrType].c_str()), "READ");
        if (!templatefile_PtDiff || !templatefile_PtDiff->IsOpen()) {
            std::cerr << "Cannot open template file: " << Form("./ProcessOutput/PtDiff/BootstrapSample_%s_%s_%d_%d_Pt_%0.1f_%0.1f_%s.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange, pTMin, pTMax, DihadronCorrTypeName[templ.corrType].c_str()) << std::endl;
            exit(1);
        }
    }

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

    VnUnit* vnResult = fitSample(isNch, templatefile, templ, datafile, data, -1);
    if (!vnResult) {
        std::cerr << "Cannot fit sample: " << data.fileNameSuffix << std::endl;
        exit(1);
    }
    VnUnit* vnResult_PtDiff = nullptr;
    if (pTMin > 0 && pTMax > 0) {
        vnResult_PtDiff = fitSample(isNch, templatefile_PtDiff, templ, datafile_PtDiff, data, -1, pTMin, pTMax);
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
        VnUnit* vnTemp = fitSample(isNch, templatefile, templ, datafile, data, sample);
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
            VnUnit* vnTemp_PtDiff = fitSample(isNch, templatefile_PtDiff, templ, datafile_PtDiff, data, sample, pTMin, pTMax);
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
VnUnit* fitSample(Bool_t isNch, TFile* templatefile, InputUnit templ, TFile* datafile, InputUnit data, int sample = -1, Double_t pTMin=0, Double_t pTMax=0) {
    std::vector<Double_t> fParamVal;
    std::vector<Double_t> fParamErr;
    TH1D* lm=0;
    TH1D* hm=0;
    TString suffix = (sample == -1) ? "" : Form("_%d", sample);
    lm = (TH1D*)templatefile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d%s", templ.minRange, templ.maxRange, suffix.Data()));
    hm = (TH1D*)datafile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d%s", data.minRange, data.maxRange, suffix.Data()));
    if (!lm) {
        std::cerr << "Cannot find histogram: " << Form("bsSample_hPhiSameOverMixed_%d_%d%s", templ.minRange, templ.maxRange, suffix.Data()) << std::endl;
        return 0;
    }
    if (!hm) {
        std::cerr << "Cannot find histogram: " << Form("bsSample_hPhiSameOverMixed_%d_%d%s", data.minRange, data.maxRange, suffix.Data()) << std::endl;
        return 0;
    }
    RooTempFitter(lm, hm, fParamVal, fParamErr, kFALSE);
    if (sample == -1) {
        PlotFitting(lm, hm, isNch, data.fileNameSuffix, data.minRange, data.maxRange, fParamVal, fParamErr, data.corrType, pTMin, pTMax);
    }
    VnUnit* vnResult = new VnUnit(fParamVal[0], fParamErr[0], fParamVal[1], fParamErr[1], fParamVal[2], fParamErr[2]);
    return vnResult;
}

//==============================================================
void RooTempFitter(TH1 *lm, TH1 *hm, std::vector<Double_t>& fParamVal, std::vector<Double_t>& fParamErr, Bool_t kRefit) {
    fParamVal.clear();
    fParamErr.clear();
    fParamVal.resize(5);
    fParamErr.resize(5);
    if (!lm ||!hm) {
        std::cerr << "Null pointer to histogram" << std::endl;
        return;
    }
    // Option: inflate or regularize per-bin errors at runtime to stabilize fits.
    // Controlled via environment variables (no code edits needed):
    // LR_ERROR_SCALE  - multiplicative factor applied to all bin errors (default 1.0)
    // LR_SYS_FRAC     - fractional systematic added in quadrature (e.g. 0.01 = 1%% of bin content)
    // LR_ERROR_FLOOR  - minimum absolute error allowed per bin (default 1e-6)
    {
        const char* es = getenv("LR_ERROR_SCALE");
        const char* sf = getenv("LR_SYS_FRAC");
        const char* ef = getenv("LR_ERROR_FLOOR");
        double envScale = es ? atof(es) : 1.0;
        double sysFrac = sf ? atof(sf) : 0.0;
        double errFloor = ef ? atof(ef) : 1e-6;
        if (!TMath::Finite(envScale) || envScale <= 0) envScale = 1.0;
        if (!TMath::Finite(sysFrac) || sysFrac < 0) sysFrac = 0.0;
        if (!TMath::Finite(errFloor) || errFloor < 0) errFloor = 1e-6;
        if (envScale != 1.0 || sysFrac != 0.0 || errFloor != 1e-6) {
            std::cout << "[FitDiag] Applying runtime error regularization: scale=" << envScale << " sysFrac=" << sysFrac << " floor=" << errFloor << std::endl;
            int nb = hm->GetNbinsX();
            for (int ib = 1; ib <= nb; ++ib) {
                double err = hm->GetBinError(ib);
                double val = hm->GetBinContent(ib);
                if (!TMath::Finite(err) || err <= 0) err = (val > 0) ? sqrt(val) : 1.0;
                if (sysFrac > 0.0) {
                    double sys = fabs(val) * sysFrac;
                    err = sqrt(err*err + sys*sys);
                }
                err *= envScale;
                if (!TMath::Finite(err) || err < errFloor) err = errFloor;
                hm->SetBinError(ib, err);
            }
        }
    }
    // Quick diagnostics: print means/max to help spot scale mismatches
    double hm_mean = hm->GetMean();
    double lm_mean = lm->GetMean();
    double hm_max = hm->GetMaximum();
    double lm_max = lm->GetMaximum();
    std::cout << "[FitDiag] pre-fit means: hm_mean=" << hm_mean << " lm_mean=" << lm_mean << " hm_max=" << hm_max << " lm_max=" << lm_max << std::endl;

    // Normalize template histogram to data mean to make Fa around O(1) where possible.
    if (TMath::Finite(lm_mean) && lm_mean != 0.0) {
        double scale_lm_to_hm = (TMath::Finite(hm_mean) && hm_mean>0) ? (hm_mean / lm_mean) : 1.0;
        if (TMath::Finite(scale_lm_to_hm) && fabs(scale_lm_to_hm - 1.0) > 1e-12) {
            std::cout << "[FitDiag] Scaling lm by " << scale_lm_to_hm << " to match hm mean" << std::endl;
            lm->Scale(scale_lm_to_hm);
        }
    }

    // Estimate Fa and Ga via a weighted linear regression: hm ~= Fa * lm + Ga
    double S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
    int nbins = hm->GetNbinsX();
    for (int ib = 1; ib <= nbins; ++ib) {
        double x = lm->GetBinContent(ib);
        double y = hm->GetBinContent(ib);
        double err = hm->GetBinError(ib);
        double w = (err > 0 && TMath::Finite(err)) ? 1.0 / (err*err) : 1.0;
        S += w;
        Sx += w * x;
        Sy += w * y;
        Sxx += w * x * x;
        Sxy += w * x * y;
    }
    double denom = (S * Sxx - Sx * Sx);
    double guessFa = 1.0;
    double guessGa = 0.0;
    if (TMath::Finite(denom) && fabs(denom) > 1e-18) {
        guessFa = (S * Sxy - Sx * Sy) / denom;
        guessGa = (Sy - guessFa * Sx) / S;
    } else {
        // fallback: use ratio of maxima or means
        if (lm_max > 0) guessFa = hm_max / lm_max;
        else if (lm_mean > 0) guessFa = hm_mean / lm_mean;
        guessGa = hm_mean - guessFa * lm_mean;
    }
    if (!TMath::Finite(guessFa)) guessFa = 1.0;
    if (!TMath::Finite(guessGa)) guessGa = 0.0;
    // Clamp guesses to reasonable ranges
    double fa_floor = 1e-6;
    if (fabs(guessFa) < fa_floor) guessFa = (guessFa >= 0) ? fa_floor : -fa_floor;
    std::cout << "[FitDiag] Regression guess: Fa=" << guessFa << " Ga=" << guessGa << std::endl;

    //Initialize fitter with given projections
    TemplateFitter *ft = new TemplateFitter(hm);
    //Setting up variable ( = delta phi, or just "x"):
    ft->AddVariable("x", "x", -TMath::Pi()/2.0, 1.5*TMath::Pi());

    if (!kRefit){
        // Pb-Pb initial value
        // OLD: ft->AddParameter("Fa","Fa",4.5,0,100);
        // OLD: ft->AddParameter("Ga","Ga",23000,0,100000);
        // Use regression-based initial guesses with data-driven bounds to improve conditioning
        {
            double fa_init = guessFa;
            // Enforce Fa non-negative: physical template scale should be >= 0.
            if (!TMath::Finite(fa_init)) fa_init = fa_floor;
            if (fa_init < 0) {
                std::cout << "[FitDiag] Regression gave negative Fa (" << guessFa << ") - using abs(fa) fallback." << std::endl;
                fa_init = fabs(fa_init);
            }
            if (fa_init < fa_floor) fa_init = fa_floor;

            double fa_min = 0.0; // do not allow negative Fa
            double fa_max = std::max(fa_min + 1.0, fa_init * 10.0);
            ft->AddParameter("Fa","Fa",fa_init,fa_min,fa_max);

            double ga_init = guessGa;
            if (!TMath::Finite(ga_init)) ga_init = 0.0;
            double ga_span = std::max(1.0, fabs(ga_init)*2.0);
            double ga_min = ga_init - ga_span;
            double ga_max = ga_init + ga_span;
            ft->AddParameter("Ga","Ga",ga_init,ga_min,ga_max);
        }
        ft->AddParameter("v2","v2",4e-3,-1.0,1.0);
        ft->AddParameter("v3","v3",6e-4,-1.0,1.0);
        ft->AddParameter("v4","v4",1.8e-4,-1.0,1.0);
        // ft->AddParameter("Fa","Fa",4.5,4.,1e3);
        // ft->AddParameter("Ga","Ga",30000,28350.,1e10);
        // ft->AddParameter("v2","v2",4e-3,0.0015,1.0);
        // ft->AddParameter("v3","v3",6e-4,0.0004,1.0);
        // ft->AddParameter("v4","v4",1.8e-4,0.00005,1.0);

        // pp initial value
        // ft->AddParameter("Fa","Fa",0.1,0,30);
        // ft->AddParameter("Ga","Ga",1,0,30);
        // ft->AddParameter("v2","v2",4e-3,-1.0,1.0);
        // ft->AddParameter("v3","v3",6e-4,-1.0,1.0);
        // ft->AddParameter("v4","v4",1.8e-4,-1.0,1.0);
        // // ft->AddParameter("v5","v5",1e-4,0,1.0);
    }
    else{
        ft->AddParameter("Fa","Fa",fParamVal[3],fParamVal[3]-1,fParamVal[3]+1);
        ft->AddParameter("Ga","Ga",fParamVal[4],fParamVal[4]-5,fParamVal[4]+5);
        ft->AddParameter("v2","v2",fParamVal[0],fParamVal[0]-0.0002,fParamVal[0]+0.0002);
        ft->AddParameter("v3","v3",fParamVal[1],fParamVal[1]-0.0002,fParamVal[1]+0.0002);
    }

    //Construct fit function
    FunctionObject *fobj = new TemplateFunction(lm);
    ft->SetFitFunction(fobj);
    //Perform fit:
    printf("About to fit\n");
    Int_t dummy = ft->Fit(0); //Do not draw performance at this point. Return value of Fit() is false if no base is set.
    if(!dummy) return;

    Double_t F   = ft->getVal(0); //0 for F, 1 for G, 2 for v2, 3 for v3, 4 for v4
    Double_t Fe  = ft->getErr(0);
    Double_t G   = ft->getVal(1); //0 for F, 1 for G, 2 for v2, 3 for v3, 4 for v4
    Double_t Ge  = ft->getErr(1);
    Double_t v21 = ft->getVal(2); //0 for F, 1 for G, 2 for v2, 3 for v3, 4 for v4
    Double_t v21e= ft->getErr(2);
    Double_t v31 = ft->getVal(3); //0 for F, 1 for G, 2 for v2, 3 for v3, 4 for v4
    Double_t v31e= ft->getErr(3);
    Double_t v41 = ft->getVal(4); //0 for F, 1 for G, 2 for v2, 3 for v3, 4 for v4
    Double_t v41e= ft->getErr(4);

    printf("Values from fit:\n");
    printf("F  = %f +- %f\n",F,Fe);
    printf("G  = %f +- %f\n",G,Ge);
    printf("V2 = %f +- %f\n",v21,v21e);
    printf("V3 = %f +- %f\n",v31,v31e);
    printf("V4 = %f +- %f\n",v41,v41e);

    fParamVal[0] = v21; fParamErr[0] = v21e;
    fParamVal[1] = v31; fParamErr[1] = v31e;
    fParamVal[2] = v41; fParamErr[2] = v41e;
    fParamVal[3] = F;   fParamErr[3] = Fe;
    fParamVal[4] = G;   fParamErr[4] = Ge;
}

void DrawText(double xmin, double ymin, double textSize, TString text)
{

	TLatex *textPreliminary = new TLatex(xmin, ymin, Form("%s", text.Data()));
	textPreliminary->SetNDC();
	textPreliminary->SetTextFont(43);
	textPreliminary->SetTextSize(textSize);
	textPreliminary->Draw();

}

void PlotFitting(TH1 *lm, TH1 *hm, Bool_t isNch, std::string fileSuffix, Int_t minRange, Int_t maxRange, const std::vector<Double_t>& par, const std::vector<Double_t>& parerr, Int_t corrType, Double_t pTMin=0, Double_t pTMax=0) {
    gStyle->SetOptStat(0); 
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    double v21 = par[0];
    double v31 = par[1];
    double v41 = par[2];
    double F =   par[3];
    double G =   par[4];

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

    // 创建拟合曲线
    const int pointBin = (int)hm->GetNbinsX();
    Double_t CopyPointX[pointBin];
    Double_t CopyPointY[pointBin];
    for (int i=0; i<pointBin; ++i){
      CopyPointX[i] = hm->GetBinCenter(i+1);
      double x = hm->GetBinCenter(i+1);
      // CopyPointY[i] = F*lm->GetBinContent(i+1)+G*(1+2*v21*cos(2*x)+2*v31*cos(3*x));
      CopyPointY[i] = F*lm->GetBinContent(i+1)+G*(1+2*v21*cos(2*x)+2*v31*cos(3*x)+2*v41*cos(4*x));
    };
    TGraph* gCopy = new TGraph(pointBin,CopyPointX,CopyPointY);

    Double_t PeriPointY[pointBin];
    for (int i=0; i<pointBin; ++i){
      PeriPointY[i] = F*lm->GetBinContent(i+1) + G;
    };
    TGraph* gPeri = new TGraph(pointBin,CopyPointX,PeriPointY);

    Double_t Y0Position = lm->GetXaxis()->FindBin(0.0);
    Double_t Y0 = lm->GetBinContent(Y0Position);


    gCopy->SetLineColor(colors[0]);
    gCopy->SetLineWidth(2);
    gCopy->Draw("same");

    gPeri->SetLineColor(colors[1]);
    gPeri->SetLineWidth(2);
    gPeri->SetLineStyle(5);
    gPeri->Draw("same");


    TF1* fit_p2 = new TF1("fit_p2","[0]*[1] + [2]*(1 + 2*[3]*cos(2*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p3 = new TF1("fit_p3","[0]*[1] + [2]*(1 + 2*[3]*cos(3*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p4 = new TF1("fit_p4","[0]*[1] + [2]*(1 + 2*[3]*cos(4*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    fit_p2->SetParameters(F,Y0,G,v21);
    fit_p3->SetParameters(F,Y0,G,v31);
    fit_p4->SetParameters(F,Y0,G,v41);
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
    leg->AddEntry(gCopy, "FY(#Delta#phi)^{peri} + G(1+#Sigma_{n=2}^{4}2V_{n#Delta}cos(n#Delta#phi))", "l");
    leg->AddEntry(gPeri, "FY(#Delta#phi)^{peri} + G", "l");
    leg->AddEntry(fit_p2, "FY(0)^{peri} + G(1+2V_{2#Delta}cos(2#Delta#phi))", "l");
    leg->AddEntry(fit_p3, "FY(0)^{peri} + G(1+2V_{3#Delta}cos(3#Delta#phi))", "l");
    leg->AddEntry(fit_p4, "FY(0)^{peri} + G(1+2V_{4#Delta}cos(4#Delta#phi))", "l");
    leg->Draw();

    // 添加文本标签
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(43);
    tex->SetTextSize(20);
    DrawText(0.2, 0.85, 20, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
    DrawText(0.2, 0.80, 20, Form("%.1f < #Delta#eta < %.1f", DihadrondEtaRange[corrType][0], DihadrondEtaRange[corrType][1]));
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

    TF1* fit_p1m = new TF1("fit_p1m","[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    fit_p1m->SetParameters(G,v21,v31,v41);

    
    // 创建残差直方图
    TH1D* hsubtract = (TH1D*)hm->Clone(Form("subtract"));
    hsubtract->Add(lm, -F);

    // =============== 新增：计算chi2/ndf ===============
    double chi2 = 0.0;
    int nBins = hsubtract->GetNbinsX();
    int nParams = 5; // 参数个数: F,G,v21,v31,v41
    int ndf = nBins - nParams;
    
    for (int i = 1; i <= nBins; i++) {
        double data = hsubtract->GetBinContent(i);
        double error = hsubtract->GetBinError(i);
        double x = hsubtract->GetBinCenter(i);
        double fit = fit_p1m->Eval(x);
        
        if (error > 0) { // 忽略误差为0的bin
            double residual = data - fit;
            chi2 += (residual * residual) / (error * error);
        }
    }
    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;
    // =============== 新增：添加chi2/ndf标签 ===============
    pad1->cd();
    TLatex* chi2Label = new TLatex();
    chi2Label->SetNDC();
    chi2Label->SetTextFont(43);
    chi2Label->SetTextSize(20);
    chi2Label->DrawLatex(0.50, 0.60, Form("#chi^{2}/ndf = %.1f/%d = %.2f", chi2, ndf, chi2ndf));

    pad2->cd();
    // TH1D* hResidual = (TH1D*)hsubtract->Clone("residual");
    // for(int i=1; i<=hResidual->GetNbinsX(); ++i){
    //   double data = hm->GetBinContent(i);
    //   double fit = fit_p1m->Eval(hm->GetBinCenter(i));
    //   hResidual->SetBinContent(i, data/fit);
    //   hResidual->SetBinError(i, hm->GetBinError(i)/fit);
    // }
        TH1D* hResidual = (TH1D*)hsubtract->Clone(Form("pull"));
        // OLD: TH1D* hResidual = (TH1D*)hsubtract->Clone(Form("pull"));
        // Guard against Inf/NaN propagated from Data/Fit division. Compute per-bin
        // Data/Fit only when the fit value is finite and non-zero. Expand the
        // display range conservatively if values are degenerate.
        double ymin_pull = 1e9, ymax_pull = -1e9;
        for (int i = 1; i <= hResidual->GetXaxis()->GetNbins(); i++) {
            double dat = hResidual->GetBinContent(i);
            double err = hResidual->GetBinError(i);
            double xcenter = hResidual->GetBinCenter(i);
            double lin = fit_p1m->Eval(xcenter);
            double ratio = 1.0;
            double ratio_err = 0.0;
            if (TMath::Finite(lin) && fabs(lin) > 1e-12) {
                ratio = dat / lin;
                if (TMath::Finite(err)) ratio_err = err / fabs(lin);
            } else {
                // If the fit value is zero/invalid, set a neutral ratio=1 and zero err
                ratio = 1.0;
                ratio_err = 0.0;
            }
            if (!TMath::Finite(ratio)) ratio = 1.0;
            if (!TMath::Finite(ratio_err)) ratio_err = 0.0;
            hResidual->SetBinContent(i, ratio);
            hResidual->SetBinError(i, ratio_err);
            if (ratio < ymin_pull) ymin_pull = ratio;
            if (ratio > ymax_pull) ymax_pull = ratio;
        }
        // If we didn't find any finite values, fall back to a safe display range
        if (!TMath::Finite(ymin_pull) || !TMath::Finite(ymax_pull) || ymax_pull <= ymin_pull) {
            ymin_pull = 0.5;
            ymax_pull = 1.5;
        } else {
            // Expand the range around 1.0 to make deviations visible, but keep a
            // minimum padding to avoid zero-height pad warnings.
            double mid = 1.0;
            double lowDelta = mid - ymin_pull;
            double highDelta = ymax_pull - mid;
            double expandLow = TMath::Max(0.1, lowDelta * 2.0);
            double expandHigh = TMath::Max(0.1, highDelta * 2.0);
            ymin_pull = mid - expandLow;
            ymax_pull = mid + expandHigh;
        }
        // Safety: ensure a non-zero delta to avoid TPad Resize warnings
        if (ymax_pull - ymin_pull < 1e-6) { ymax_pull = ymin_pull + 1e-3; }
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
        canvas->SaveAs(Form("./TemplateFit/PtDiff/PDFs/TemplateFit_%s_%s_%d_%d_Pt_%0.1f_%0.1f_%s.pdf", fileSuffix.c_str(), splitName.c_str(), minRange, maxRange, pTMin, pTMax, DihadronCorrTypeName[corrType].c_str()));
    } else {
        canvas->SaveAs(Form("./TemplateFit/PDFs/TemplateFit_%s_%s_%d_%d_%s.root", fileSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
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

