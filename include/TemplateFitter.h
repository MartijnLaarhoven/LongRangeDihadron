/*
Vytautas Vislavicius
vytautas.vislavicius@cern.ch

Generic 1/2/3D chi2 fitter using RooFit. Functions, templates, you name it, it fits it.
Code is simple, but a lot of thought went into it. I have no way to control who, how,
and for which purpose one uses it, so please be fair and whenever using it, at least
acknowledge clearly the author of this.
*/
#ifndef TEMPLATEFITTER__H
#define TEMPLATEFITTER__H
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"
#include "RooRealSumPdf.h"
#include "RooParamHistFunc.h"
#include "RooHistConstraint.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TNamed.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "RooTFnBinding.h"
#include "FunctionObject.h"
#include "TRandom.h"
#include "TProfile.h"
using namespace RooFit;
class TemplateFitter {
 public:
    TemplateFitter();
    TemplateFitter(TH1 *datah): dataH(0),dataBU(0),fStatErr(0),totFunc(0),fParList(0),fVarList(0),f_Dim(0),f_FObj(0),f_FitFunc(0),f_Ready(0) { SetData(datah); };
    ~TemplateFitter();
    void SetData(TH1 *inh);
    void Reset() { if(dataH) printf("Warning! This will delete the target histogram! Call SetData(TH1*) afterwards!\n"); delete dataH; f_Ready=kFALSE; };
    void SetFitFunction(FunctionObject *fobj);
    void AddParameter(TString lName, TString lTitle, Double_t l_val, Double_t l_min, Double_t l_max); //parameters
    void AddVariable(TString lName, TString lTitle, Double_t l_min, Double_t l_max); //variables ( = dimensions)
    Bool_t SetConst(Int_t ind, Bool_t toConst);
    Bool_t SetConst(Int_t ind, Double_t cVal);
    void Clear() { fParList->Clear(); fVarList->Clear(); delete totFunc; totFunc=0; delete dataH; dataH=0;};
    Bool_t PrepareForFitting();
    TList *getParList() { return fParList; };
    // void ResetParList() {if(fArgList) delete fArgList; };
    // void SetBase(TH1 *inh) {if(tmplH) delete tmplH; tmplH = (TH2*)inh->Clone("l_tmplH"); tmplH->SetDirectory(0); };
    Bool_t Fit(Int_t nRefits=0);
    Double_t getVal(Int_t ind) {return ((RooRealVar*)fParList->At(ind))->getVal();}; //0 for G, 1 for F, 2 for v2, 3 for v3
    Double_t getErr(Int_t ind) {return ((RooRealVar*)fParList->At(ind))->getError();};
    TF1 *getFitFunction() { return totFunc; };
    void Randomize(Bool_t RestoreOriginal=kFALSE);
    TProfile *getStatErrors() { if(!fStatErr->GetEntries()) printf("Statistics are not profiled! To profile, run Fit(N>0), where N is a number of itterations\n"); return fStatErr; };
  protected:
    TH1 *dataH;
    TH1 *dataBU;
    TProfile *fStatErr;
    TF1 *totFunc;
    TList *fParList;
    TList *fVarList;
    Int_t f_Dim;
    FunctionObject *f_FObj; //Shallow copy
    RooAbsReal *f_FitFunc;
    Bool_t f_Ready;
    Int_t getDimension();
    Bool_t rescaleHistogram(Bool_t divide=kFALSE);
    Int_t getAxNbins(TAxis *inax) {if(!inax) return 0; return inax->GetXbins()->fN; };
    void ScaleBin(TH1 *inh, Double_t scale, Int_t bx, Int_t by=0, Int_t bz=0);
    Bool_t SetupFF();
};
#endif
