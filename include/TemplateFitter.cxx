/*
Vytautas Vislavicius
vytautas.vislavicius@cern.ch

Generic 1/2/3D chi2 fitter using RooFit. Functions, templates, you name it, it fits it.
Code is simple, but a lot of thought went into it. I have no way to control who, how,
and for which purpose one uses it, so please be fair and whenever using it, at least
acknowledge clearly the author of this.
*/
#include "TemplateFitter.h"
// #include "RooChi2Var.h"
using namespace RooFit;
TemplateFitter::TemplateFitter():
  dataH(0),
  dataBU(0),
  totFunc(0),
  fParList(0)
{
};
TemplateFitter::~TemplateFitter() {
  delete dataH;
  delete dataBU;
  delete totFunc;
  delete fParList;
};
void TemplateFitter::SetData(TH1 *inh) {
   if(dataH) delete dataH;
   dataH = (TH1*)inh->Clone("l_dataH");
   dataH->SetDirectory(0);
   f_Dim=getDimension();
   f_Ready=kFALSE;
   dataBU = (TH1*)dataH->Clone("l_dataH_Backup");
   dataBU->SetDirectory(0);
}
Bool_t TemplateFitter::PrepareForFitting() {
  if(f_Ready) return kTRUE;
  if(!f_Dim) {printf("Input histogram has dimension 0, cannot fit...\n"); return 0; };
  if(!fParList || !fParList->GetEntries()) {printf("No parameters specified. You probably forgot to call AddParameter().\n"); return 0; };
  Int_t l_nDim = fVarList->GetEntries();
  if(l_nDim!=f_Dim) {printf("Number of variables (%i) is not the same as number of histogram dimensions (%i). Won't fit.\n",l_nDim,f_Dim); return 0;};
  if(!rescaleHistogram(kFALSE)) {printf("Could not rescale the histogram. Quitting...\n"); return 0; };
  //Figure out dimensions & define variables
  RooRealVar *x, *y, *z;
  Double_t xmin=0,xmax=0,ymin=0,ymax=0,zmin=0,zmax=0;
  if(l_nDim>0) {
    x = ((RooRealVar*)fVarList->At(0));
    xmin=x->getMin(); xmax=x->getMax();
  }
  if(l_nDim>1) {
    y = ((RooRealVar*)fVarList->At(1));
    ymin=y->getMin(); ymax=y->getMax();
  };
  if(l_nDim>2) {
    z = ((RooRealVar*)fVarList->At(2));
    zmin=z->getMin(); zmax=z->getMax();
  };
  if(!SetupFF()) return kFALSE;
  // totFunc->SetRange(xmin,xmax,ymin,ymax,zmin,zmax);
  totFunc->SetRange(xmin,xmax);
  RooArgList lArgList;
  for(Int_t i=0;i<fParList->GetEntries(); i++) lArgList.add(*((RooRealVar*)fParList->At(i)));
  //Set up the function
  if(f_FitFunc) delete f_FitFunc;
  if(l_nDim==1) f_FitFunc = bindFunction((TF1*)totFunc,*x,lArgList);
  else if(l_nDim==2) f_FitFunc = bindFunction((TF2*)totFunc,*x,*y,lArgList);
  else if(l_nDim==3) f_FitFunc = bindFunction((TF3*)totFunc,*x,*y,*z,lArgList);
  else { printf("Currently, only 1-3 dimensions are supported.\n"); return kFALSE; };
  //Prepare profiling for statistics
  if(!fStatErr) delete fStatErr;
  fStatErr = new TProfile("StatErrors","Stat errors",fParList->GetEntries(),0,fParList->GetEntries());
  fStatErr->SetErrorOption("s");
  for(Int_t i=0; i<fParList->GetEntries(); i++)
    fStatErr->GetXaxis()->SetBinLabel(i+1,((RooRealVar*)fParList->At(i))->GetName());
  f_Ready = kTRUE;
  return kTRUE;
}
Bool_t TemplateFitter::Fit(Int_t nRefits) {
  if(!PrepareForFitting()) return 0;
  Int_t l_nDim = fVarList->GetEntries();
  RooArgList lVarList;
  for(Int_t i=0;i<l_nDim;i++) lVarList.add(*((RooRealVar*)fVarList->At(i)));
  RooDataHist dsig("dsig", "dsig", lVarList, Import(*dataH,kFALSE));
  //Perform fit
  f_FitFunc->chi2FitTo(dsig, RooFit::PrintLevel(-1), RooFit::SumW2Error(kTRUE), RooFit::Warnings(kFALSE));//signal is divided by bin width at this point


  // Check chi2
  // RooAbsReal* chi2_o_ndf = f_FitFunc->createChi2(dsig, Range("cutrange"), Extended(true), DataError(RooAbsData::Poisson));
  RooAbsReal* chi2_o_ndf = f_FitFunc->createChi2(dsig, Range("cutrange"), Extended(true), DataError(RooAbsData::SumW2));
  std::cout << "chi2: " << chi2_o_ndf->getVal() << std::endl;
  const int ndf = (int)dataH->GetNbinsX()-5;
  float chi2ndf = chi2_o_ndf->getVal() / ndf;
  cout<< "ndf: " << ndf << ", chi2/ndf show  :  " << chi2ndf <<endl;
  // reject if chi2/ndf is too large
  if (chi2ndf > 30) {
    // cout << "Chi2/ndf is too large, setting all parameters to -1" << endl;
    for(Int_t i=0;i<fParList->GetEntries();i++) {
      ((RooRealVar*)fParList->At(i))->setVal(-1);
      ((RooRealVar*)fParList->At(i))->setError(10);
    };
    return 1;
  }


  // // Draw and check fit // ediited by Wenya, Jan 2023
  // cout<<"now trying to draw fit quality plots!"<<endl;
  // TCanvas* c2 = new TCanvas();
  // RooRealVar *x;
  // Double_t xmin=0,xmax=0;
  // if(l_nDim>0) {
  //   x = ((RooRealVar*)fVarList->At(0));
  //   xmin=x->getMin(); xmax=x->getMax();
  // }
  // RooDataHist dsigBU("dsigBU", "dsigBU", lVarList, Import(*dataBU,kFALSE));
  // RooPlot* frame = x->frame();
  // dsigBU.plotOn(frame, LineColor(kRed), MarkerColor(kRed));
  // f_FitFunc->plotOn(frame,  MarkerColor(kBlue));
  // c2->cd();
  // frame->Draw();
  // totFunc->Draw("same l");
  // c2->SaveAs("/home/liscklu/Workdir/O2Tutorial/DihardonCorr/temp/fit_quality.pdf");
  // cout<<"End draw!"<<endl;


  if(nRefits<=0) return 1;



  // for(Int_t i=0;i<fParList->GetEntries();i++) {mws.push_back(0); ws.push_back(0);};
  for(Int_t i=0;i<nRefits;i++) {
    f_FObj->Randomize();
    Randomize(kFALSE);
    f_FitFunc->chi2FitTo(dsig, RooFit::PrintLevel(-1), RooFit::SumW2Error(kTRUE), RooFit::Warnings(kFALSE));
    for(Int_t j=0;j<fParList->GetEntries();j++)
      fStatErr->Fill(j,getVal(j),getErr(j));
  };

  for(Int_t i=0;i<fParList->GetEntries();i++) {
    ((RooRealVar*)fParList->At(i))->setVal(fStatErr->GetBinContent(i+1));
    ((RooRealVar*)fParList->At(i))->setError(fStatErr->GetBinError(i+1));
  };




  Randomize(kTRUE); //restore the original data histogram
  f_FObj->Restore();

  return 1;
};
void TemplateFitter::SetFitFunction(FunctionObject *fobj) {
  if(!fobj) {printf("Null pointer for function object passed! Not doing anything this time...\n"); return; };
  if(!fobj->isValid()) {printf("Fit function is not valid! Check isValid() method in your function implementation...\n"); return; };
  f_FObj = fobj; //Just a shallow copy
};
void TemplateFitter::AddParameter(TString lName, TString lTitle, Double_t l_val, Double_t l_min, Double_t l_max) {
  if(!fParList) {
    fParList = new TList();
    fParList->SetOwner(kTRUE);
  };
  RooRealVar *l_Var = new RooRealVar(lName,lTitle,l_val,l_min,l_max);
  fParList->Add(l_Var);
}
void TemplateFitter::AddVariable(TString lName, TString lTitle, Double_t l_min, Double_t l_max) {
  if(!fVarList) {
    fVarList = new TList();
    fVarList->SetOwner(kTRUE);
  };
  RooRealVar *l_Var = new RooRealVar(lName,lTitle,l_min,l_max);
  fVarList->Add(l_Var);
}
Bool_t TemplateFitter::SetConst(Int_t ind, Bool_t isConst) {
  if(ind>=fParList->GetEntries()) {
    printf("Could not find the parameter requested (asked for var %i, while there are only %i variables)!\n",ind,fParList->GetEntries());
    if(ind==fParList->GetEntries()) printf("Are you sure you did not forget that parameters start with index 0?\n");
    return kFALSE;
  };
  RooRealVar *vr = (RooRealVar*)fParList->At(ind);
  vr->setConstant(isConst);
  return kTRUE;
};
Bool_t TemplateFitter::SetConst(Int_t ind, Double_t cVal) {
  if(ind>=fParList->GetEntries()) {
    printf("Could not find the parameter requested (asked for var %i, while there are only %i variables)!\n",ind,fParList->GetEntries());
    if(ind==fParList->GetEntries()) printf("Are you sure you did not forget that parameters start with index 0?\n");
    return kFALSE;
  };
  RooRealVar *vr = (RooRealVar*)fParList->At(ind);
  vr->setVal(cVal);
  vr->setConstant(kTRUE);
  return kTRUE;
};
Int_t TemplateFitter::getDimension() {
  if(!dataH) return 0;
  Int_t rDim = 0;
  if(dynamic_cast<TH3*>(dataH)) rDim=3;
  else if(dynamic_cast<TH2*>(dataH)) rDim=2;
  else if(dynamic_cast<TH1*>(dataH)) rDim=1;
  return rDim;
}
Bool_t TemplateFitter::rescaleHistogram(Bool_t divide) {
  if(!f_Dim) getDimension();
  if(!f_Dim) return 0;
  Int_t fNx = getAxNbins(dataH->GetXaxis());
  Int_t fNy = getAxNbins(dataH->GetYaxis());
  Int_t fNz = getAxNbins(dataH->GetZaxis());
  Double_t scale = 1;
  if(!fNx) scale*=dataH->GetXaxis()->GetBinWidth(1);
  if(!fNy && f_Dim>1) scale*=dataH->GetYaxis()->GetBinWidth(1);
  if(!fNz && f_Dim>2) scale*=dataH->GetZaxis()->GetBinWidth(1);
  if(!fNx && !fNy && !fNz) { dataH->Scale(divide?(1./scale):scale); return kTRUE; }; //In case all bin are homogenious in their dimensions, then just scale the whole histogram
  for(Int_t i=1;i<=dataH->GetNbinsX();i++) {
    Double_t x_scale=dataH->GetXaxis()->GetBinWidth(i);
    if(f_Dim<2) ScaleBin(dataH,divide?(1./x_scale):x_scale,i);
    else {
      for(Int_t j=1;j<=dataH->GetNbinsY();j++){
        Double_t xy_scale=x_scale*dataH->GetYaxis()->GetBinWidth(j);
        if(f_Dim<3) ScaleBin(dataH,divide?(1./xy_scale):xy_scale,i,j);
        else {
          for(Int_t k=1;k<=dataH->GetNbinsZ();k++) {
            Double_t xyz_scale=xy_scale*dataH->GetZaxis()->GetBinWidth(k);
            ScaleBin(dataH,divide?(1./xyz_scale):xyz_scale,i,j,k);
          }
        }
      }
    }
  }
  return kTRUE;
}
void TemplateFitter::ScaleBin(TH1 *inh, Double_t scale, Int_t bx, Int_t by, Int_t bz) {
  Int_t binNo = inh->GetBin(bx,by,bz);
  inh->SetBinContent(binNo,inh->GetBinContent(binNo)*scale);
  inh->SetBinError(binNo,inh->GetBinError(binNo)*scale);
}
Bool_t TemplateFitter::SetupFF() {
  if(!f_FObj) {printf("FunctionObject not set! You probably forgot to call SetFitFunction()...\n"); return kFALSE; };
  if(totFunc) delete totFunc;
  totFunc = new TF1("total_FitFunction",f_FObj,0,10,fParList->GetEntries());
  return kTRUE;
}
void TemplateFitter::Randomize(Bool_t RestoreOriginal) {
  if(!f_Ready) {
    printf("Warning! The fitter has not been called yet, and so the data histogram has not been scaled yet. Please run TemplateFitter::Fit() or TemplateFitter::PrepareForFitting() to make sure that the histogram is ready for fitting!\n");
  };
  if(!dataBU) {
    dataBU = (TH1*)dataH->Clone("l_dataH_Backup");
    dataBU->SetDirectory(0);
    if(RestoreOriginal) return;
  };
  if(RestoreOriginal) {
    if(dataH) delete dataH;
    dataH = (TH1*)dataBU->Clone("l_dataH");
    dataH->SetDirectory(0);
    return;
  }
  for(Int_t i=1;i<=dataH->GetNbinsX();i++) {
    dataH->SetBinContent(i,gRandom->Gaus(dataBU->GetBinContent(i),dataBU->GetBinError(i)));
  };
}
