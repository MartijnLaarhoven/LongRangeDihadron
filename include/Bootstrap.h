#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

#include <vector>
#include <TMath.h>
#include <iostream>

double epsilon = 1e-10;

bool isZero(double x) {
  return (x < epsilon && x > -epsilon);
}

void CalculateBootstrapError(const std::vector<std::vector<double>>& ValueArray, const std::vector<std::vector<double>>& ValueErrorArray, std::vector<double>& ErrorArray, Double_t TotalOverSubSample){
    int Nsample = ValueArray.size();
    int Nbin = ValueArray[0].size();
    std::vector<int> Count;
    std::vector<double> Mean;
    std::vector<double> SumWeight;
    std::vector<std::vector<bool>> MaskArray;
    Count.resize(Nbin);
    Mean.resize(Nbin);
    SumWeight.resize(Nbin);
    MaskArray.resize(Nsample);
    // first calculate the mean and standard deviation for each bin
    for(int isample=0;isample<Nsample;isample++){
        MaskArray[isample].resize(Nbin);
        for(int ibin=0;ibin<Nbin;ibin++){
            if (ValueArray[isample][ibin]==-1 && ValueErrorArray[isample][ibin]==10) {
                MaskArray[isample][ibin] = false;
                continue;
            }
            if (isZero(ValueArray[isample][ibin]) || isZero(ValueErrorArray[isample][ibin])) {
                MaskArray[isample][ibin] = false;
                continue;
            }
            if (ValueArray[isample][ibin]!=ValueArray[isample][ibin] || ValueErrorArray[isample][ibin]!=ValueErrorArray[isample][ibin]) {
                MaskArray[isample][ibin] = false;
                continue;
            }
            MaskArray[isample][ibin] = true;
        }
    }
    for(int ibin=0;ibin<Nbin;ibin++){
        SumWeight[ibin]=0;
        Mean[ibin]=0;
        for(int isample=0;isample<Nsample;isample++){
            if (!MaskArray[isample][ibin]) continue;
            Count[ibin]++;
            Mean[ibin] += ValueArray[isample][ibin] * (1./(ValueErrorArray[isample][ibin]*ValueErrorArray[isample][ibin]));
            // Printf("ValueArray[%d][%d]=%f, ValueErrorArray[%d][%d]=%f",j,i,ValueArray[j][i],j,i,ValueErrorArray[j][i]);
            SumWeight[ibin] += (1./(ValueErrorArray[isample][ibin]*ValueErrorArray[isample][ibin]));
        }
        if(Count[ibin]>0)
            Mean[ibin] = Mean[ibin]/SumWeight[ibin];
    }
    for (Int_t isample = 0; isample < Nsample; isample++)
    {
        for (Int_t ibin = 0; ibin < Nbin; ibin++)
        {
            if (!MaskArray[isample][ibin]) continue;
            ErrorArray[ibin] += (ValueArray[isample][ibin] - Mean[ibin]) * (ValueArray[isample][ibin] - Mean[ibin]);
        }
    }
    for (Int_t ibin = 0; ibin < Nbin; ibin++)
    {
        if (Count[ibin] > 2){
            ErrorArray[ibin] = TMath::Sqrt(ErrorArray[ibin] / (Count[ibin] - 1));
        }
        else {
            ErrorArray[ibin] = 10;
            Printf("Bin %d has only %d samples, set error to 10",ibin,Count[ibin]);
        }
        // Printf("ErrorArray[%d]=%f",i,ErrorArray[i]);
    }

    // remove samples beyond 3 sigma
    for(int ibin=0;ibin<Nbin;ibin++) {
        for(int isample=0;isample<Nsample;isample++) {
            if (!MaskArray[isample][ibin]) continue;
            if (TMath::Abs(ValueArray[isample][ibin] - Mean[ibin]) > 3 * ErrorArray[ibin]) {
                MaskArray[isample][ibin] = false;
            }
        }
    }

    // re-calculate the mean and standard deviation for each bin
    for(int ibin=0;ibin<Nbin;ibin++){
        SumWeight[ibin]=0;
        Mean[ibin]=0;
        Count[ibin]=0;
        ErrorArray[ibin]=0;
        for(int isample=0;isample<Nsample;isample++){
            if (!MaskArray[isample][ibin]) continue;
            Count[ibin]++;
            Mean[ibin] += ValueArray[isample][ibin] * (1./(ValueErrorArray[isample][ibin]*ValueErrorArray[isample][ibin]));
            SumWeight[ibin] += (1./(ValueErrorArray[isample][ibin]*ValueErrorArray[isample][ibin]));
        }
        if(Count[ibin]>0)
            Mean[ibin] = Mean[ibin]/SumWeight[ibin];
    }
    for (Int_t isample = 0; isample < Nsample; isample++)
    {
        for (Int_t ibin = 0; ibin < Nbin; ibin++)
        {
            if (!MaskArray[isample][ibin]) continue;
            ErrorArray[ibin] += (ValueArray[isample][ibin] - Mean[ibin]) * (ValueArray[isample][ibin] - Mean[ibin]);
        }
    }
    for (Int_t ibin = 0; ibin < Nbin; ibin++)
    {
        if (Count[ibin] > 2){
            ErrorArray[ibin] = TMath::Sqrt(ErrorArray[ibin] / (Count[ibin] - 1));
        }
        else {
            ErrorArray[ibin] = 10;
            Printf("Bin %d has only %d samples, set error to 10",ibin,Count[ibin]);
        }
        // Printf("ErrorArray[%d]=%f",i,ErrorArray[i]);
    }
    // SEM = SD/sqrt(N)
    // because ErrorArray[ibin] above is the error of subsample, which is smaller than the total sample,
    for (Int_t ibin = 0; ibin < Nbin; ibin++)
    {
        ErrorArray[ibin] = ErrorArray[ibin] / TMath::Sqrt(TotalOverSubSample);
    }
    
}

void ResizeValueArray(std::vector<std::vector<std::vector<double>>>& ValueArray,
    std::vector<std::vector<std::vector<double>>>& ValueErrorArray,
    std::vector<std::vector<double>>& ErrorArray,
    int Nobs=1, int NofSample=10, int Nbin=9){
    ValueArray.clear();
    ValueErrorArray.clear();
    ErrorArray.clear();
    ValueArray.resize(Nobs);
    ValueErrorArray.resize(Nobs);
    ErrorArray.resize(Nobs);    
    for(int i=0;i<Nobs;i++){
        ErrorArray[i].resize(Nbin);
        ValueArray[i].resize(NofSample);
        ValueErrorArray[i].resize(NofSample);
        for(int j=0;j<NofSample;j++)ValueArray[i][j].resize(Nbin);
        for(int j=0;j<NofSample;j++)ValueErrorArray[i][j].resize(Nbin);
    }

}


#endif // BOOTSTRAP_H
