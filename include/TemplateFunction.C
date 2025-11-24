#ifndef TEMPLATE_FUNCTION__C
#define TEMPLATE_FUNCTION__C
#include "FunctionObject.h"
#include "TMath.h"
#include "TH1.h"
#include "TRandom.h"
class TemplateFunction : public FunctionObject
{
public:
  TemplateFunction() : baseTemplate(0), backupTemplate(0) {};
  TemplateFunction(TH1 *inh) : baseTemplate(inh) {};
  ~TemplateFunction() { delete backupTemplate; };
  double operator()(double *x, double *p);
  double getError(double *x, double *p);
  Bool_t isValid() { return baseTemplate ? kTRUE : kFALSE; };
  void SetBaseTemplate(TH1 *inh) { baseTemplate = inh; };
  void Randomize();
  void Restore();

protected:
  TH1 *baseTemplate; //! Not stored by the function, do not delete
  TH1 *backupTemplate;
};
double TemplateFunction::operator()(double *x, double *p)
{
  Double_t bc = baseTemplate->GetBinContent(baseTemplate->FindBin(x[0]));
  // return p[0]*bc + p[1]*(1 + 2*p[2]*TMath::Cos(2*x[0]) + 2*p[3]*TMath::Cos(3*x[0]) );
  return p[0] * bc + p[1] * (1 + 2 * p[2] * TMath::Cos(2 * x[0]) + 2 * p[3] * TMath::Cos(3 * x[0]) + 2 * p[4] * TMath::Cos(4 * x[0]));
  // return p[0]*bc + p[1]*(1 + 2*p[2]*TMath::Cos(2*x[0]) + 2*p[3]*TMath::Cos(3*x[0]) + 2*p[4]*TMath::Cos(4*x[0]) + 2*p[5]*TMath::Cos(6*x[0]));
  // return p[0] * bc + p[1] * (1 + 2 * p[5] * TMath::Cos(x[0]) + 2 * p[2] * TMath::Cos(2 * x[0]) + 2 * p[3] * TMath::Cos(3 * x[0]) + 2 * p[4] * TMath::Cos(4 * x[0]));
};
double TemplateFunction::getError(double *x, double *p)
{ // This is absolutely optional, implemented here b/c I used the same function for TMinuit & TFractionFitters
  return baseTemplate->GetBinError(baseTemplate->FindBin(x[0])) * p[0];
};
void TemplateFunction::Randomize()
{
  backupTemplate = (TH1 *)baseTemplate->Clone(Form("%s_backup", baseTemplate->GetName()));
  for (Int_t i = 1; i <= baseTemplate->GetNbinsX(); i++)
    baseTemplate->SetBinContent(i, gRandom->Gaus(backupTemplate->GetBinContent(i), backupTemplate->GetBinError(i)));
};
void TemplateFunction::Restore()
{
  if (!backupTemplate)
  {
    printf("Nothing to restore\n");
    return;
  };
  for (Int_t i = 1; i <= baseTemplate->GetNbinsX(); i++)
  {
    baseTemplate->SetBinContent(i, backupTemplate->GetBinContent(i));
    baseTemplate->SetBinError(i, backupTemplate->GetBinError(i));
  };
}
#endif
