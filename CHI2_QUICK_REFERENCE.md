# Chi² Calculation - Quick Reference Guide

## Quick File Navigation

### Main Analysis Files

| File | Purpose | Has Chi²? | Has Syst? |
|------|---------|-----------|-----------|
| [Process_dPhidEta.cxx](Process_dPhidEta.cxx) | Extract correlations from ROOT files | ❌ | ❌ |
| [Process_CreateBootstrapSample.cxx](Process_CreateBootstrapSample.cxx) | Create bootstrap samples | ❌ | ❌ |
| [Process_FourierFit.cxx](Process_FourierFit.cxx) | Fourier fit method | ✅ | ❌ |
| [Process_TemplateFit.cxx](Process_TemplateFit.cxx) | Template fit method | ✅ | ✅ |
| [Process_3times2PC.cxx](Process_3times2PC.cxx) | 3×2PC method | ❌ | ❌ |

### Support Files

| File | Purpose |
|------|---------|
| [include/TemplateFitter.cxx](include/TemplateFitter.cxx) | Generic RooFit chi² fitter |
| [include/ErrorPropagation.h](include/ErrorPropagation.h) | Error propagation formulas |
| [include/Bootstrap.h](include/Bootstrap.h) | Bootstrap error calculation |
| [include/BasicForDihadron.h](include/BasicForDihadron.h) | Common utilities |

---

## Chi² Calculation Code Locations

### 1. Fourier Fit - RooFit Native Chi² 

**File**: [Process_FourierFit.cxx](Process_FourierFit.cxx)  
**Function**: `void RooFourierFit(TH1 *hist, ...)`  
**Start Line**: ~655  
**Chi² Calculation**: Lines 681-696

```cpp
// Line ~681
// 计算拟合质量指标
double chi2 = fitFunc->GetChisquare();
double ndf = fitFunc->GetNDF();
double chi2_ndf = (ndf > 0) ? chi2 / ndf : 0;
std::cout << "Chi2/NDF = " << chi2_ndf << std::endl;
double pvalue = TMath::Prob(chi2, static_cast<int>(ndf));
std::cout << "P-value = " << pvalue << std::endl;

if (chi2_ndf > 2000) {
    std::cout << "WARNING: Chi2/NDF > 20, fit may be inaccurate!" << std::endl;
    // Set parameters to -1 (failed fit)
    for (int i = 0; i < 5; ++i) {
        fParamVal[i] = -1;
        fParamErr[i] = 10;
    }
    return;
}
```

**Key Points**:
- ❌ NO systematic uncertainties
- Uses ROOT's native `TF1::GetChisquare()`
- Threshold: χ²/ndf > 2000 → REJECT
- Fit function: f(φ) = a₀ + 2Σaₙcos(nφ)

---

### 2. Fourier Fit - Manual Chi² for Plotting

**File**: [Process_FourierFit.cxx](Process_FourierFit.cxx)  
**Function**: `void PlotFitting(TH1 *hm, ...)`  
**Start Line**: ~833  
**Chi² Calculation**: Lines 838-855

```cpp
// Line ~838
// =============== 新增：计算chi2/ndf ===============
double chi2 = 0.0;
int nBins = hm->GetNbinsX();
int nParams = 5; // a0, a1, a2, a3, a4
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
```

**Display on Plot**: Line ~925
```cpp
chi2Label->DrawLatex(0.50, 0.60, 
    Form("#chi^{2}/ndf = %.1f/%d = %.2f", chi2, ndf, chi2ndf));
```

---

### 3. Template Fit - Systematic Uncertainty Setup

**File**: [Process_TemplateFit.cxx](Process_TemplateFit.cxx)  
**Function**: `void RooTemplateFit(TH1 *lm, TH1 *hm, ...)`  
**Start Line**: ~728  
**Systematic Setup**: Lines 748-776

```cpp
// Line ~748
// Option: inflate or regularize per-bin errors at runtime to stabilize fits.
// Controlled via environment variables (no code edits needed):
// LR_ERROR_SCALE  - multiplicative factor applied to all bin errors (default 1.0)
// LR_SYS_FRAC     - fractional systematic added in quadrature (e.g. 0.01 = 1%)
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
        std::cout << "[FitDiag] Applying runtime error regularization: "
                  << "scale=" << envScale 
                  << " sysFrac=" << sysFrac 
                  << " floor=" << errFloor << std::endl;
        
        int nb = hm->GetNbinsX();
        for (int ib = 1; ib <= nb; ++ib) {
            double err = hm->GetBinError(ib);
            double val = hm->GetBinContent(ib);
            if (!TMath::Finite(err) || err <= 0) 
                err = (val > 0) ? sqrt(val) : 1.0;
            
            // Add systematic in quadrature
            if (sysFrac > 0.0) {
                double sys = fabs(val) * sysFrac;
                err = sqrt(err*err + sys*sys);
            }
            
            err *= envScale;
            if (!TMath::Finite(err) || err < errFloor) 
                err = errFloor;
            
            hm->SetBinError(ib, err);
        }
    }
}
```

**Formula**:
```
σ_total = max(LR_ERROR_SCALE × √(σ_stat² + (LR_SYS_FRAC × |value|)²), LR_ERROR_FLOOR)
```

---

### 4. Template Fit - RooFit Chi² Calculation

**File**: [Process_TemplateFit.cxx](Process_TemplateFit.cxx)  
**Function**: `void RooTemplateFit(TH1 *lm, TH1 *hm, ...)`  
**Chi² Calculation**: Lines ~970-990 (within TemplateFitter usage)

The actual chi² is calculated inside the TemplateFitter class:

**File**: [include/TemplateFitter.cxx](include/TemplateFitter.cxx)  
**Function**: `Bool_t TemplateFitter::Fit(Int_t nRefits)`  
**Lines**: 80-110

```cpp
// Line ~80
// Perform fit
f_FitFunc->chi2FitTo(dsig, 
                     RooFit::PrintLevel(-1), 
                     RooFit::SumW2Error(kTRUE), 
                     RooFit::Warnings(kFALSE));

// Check chi2
RooAbsReal* chi2_o_ndf = f_FitFunc->createChi2(dsig, 
                                                Range("cutrange"), 
                                                Extended(true), 
                                                DataError(RooAbsData::SumW2));
std::cout << "chi2: " << chi2_o_ndf->getVal() << std::endl;
const int ndf = (int)dataH->GetNbinsX()-5;
float chi2ndf = chi2_o_ndf->getVal() / ndf;
cout<< "ndf: " << ndf << ", chi2/ndf show  :  " << chi2ndf <<endl;

// reject if chi2/ndf is too large
if (chi2ndf > 30000) {
    // cout << "Chi2/ndf is too large, setting all parameters to -1" << endl;
    for(Int_t i=0;i<fParList->GetEntries();i++) {
        ((RooRealVar*)fParList->At(i))->setVal(-1);
        ((RooRealVar*)fParList->At(i))->setError(10);
    };
    return 1;
}
```

**Key Points**:
- ✅ Uses errors modified by systematic uncertainty setup
- Uses RooFit's `createChi2()` with SumW2 error model
- Threshold: χ²/ndf > 30000 → REJECT

---

### 5. Template Fit - Manual Chi² for Plotting

**File**: [Process_TemplateFit.cxx](Process_TemplateFit.cxx)  
**Function**: `void PlotFitting(TH1 *lm, TH1 *hm, ...)`  
**Start Line**: ~1075  
**Chi² Calculation**: Lines 1088-1105

```cpp
// Line ~1085
TH1D* hsubtract = (TH1D*)hm->Clone(Form("subtract"));
hsubtract->Add(lm, -F);  // Subtract template

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
```

**Display on Plot**: Line ~1112
```cpp
chi2Label->DrawLatex(0.50, 0.60, 
    Form("#chi^{2}/ndf = %.1f/%d = %.2f", chi2, ndf, chi2ndf));
```

**Key Difference**: Fits the **background-subtracted** histogram (hsubtract), not raw data.

---

## Error Propagation Functions

**File**: [include/ErrorPropagation.h](include/ErrorPropagation.h)

### Function Reference Table

| Line | Function | Formula | Use Case |
|------|----------|---------|----------|
| 10-26 | `Error_Ratio()` | z = x/y | General ratios |
| 28-36 | `Error_vNL()` | v = N/√D | V_n from correlations |
| 38-46 | `Error_Rho()` | ρ = N/√(D₁D₂) | Correlation coefficient |
| 48-54 | `Error_Chi()` | χ = N/D | Simple ratio (wrapper) |
| 56-64 | `Error_SCnm()` | SC = N - N₁N₂ | Cumulant subtraction |
| 66-84 | `Error_SCklm()` | 3-particle correlations | Multi-particle |
| 184-192 | `Error_Ratio_sqrtY()` | z = x/√y | V_n PtDiff |

---

### Most Important Functions (with Code)

#### 1. Error_Ratio - General Ratio Error

```cpp
// Line 10-26
double Error_Ratio(double x, double ex, double y, double ey, double rho){
    //z=x/y
    //rho is correlation term, which is sqrt(N_y/N_x)
    //Cov(x,y) = rho*ex*ey
    double Contain = pow(ex/y,2)
        +pow(x*ey/(y*y),2)
        +2*(1./y)*(-1.*x/(y*y))*rho*ex*ey;
    if(Contain<0){
        rho=1;
        Contain = pow(ex/y,2)
        +pow(x*ey/(y*y),2)
        +2*(1./y)*(-1.*x/(y*y))*rho*ex*ey;
        if(Contain<0)Contain=4;
    }
    return sqrt(Contain);
}
```

**Mathematical Formula**:
$$\sigma_z = \sqrt{\left(\frac{\sigma_x}{y}\right)^2 + \left(\frac{x\sigma_y}{y^2}\right)^2 + 2\rho\frac{\sigma_x\sigma_y}{y^2}}$$

**Used for**: Converting v²ₙ to vₙ, calculating vₙ²/a₀ ratios

---

#### 2. Error_Ratio_sqrtY - Ratio with Square Root

```cpp
// Line 184-192
double Error_Ratio_sqrtY(double x, double ex, double y, double ey) {
    //z=x/sqrt(y)
    double Contain = pow(ex/sqrt(y),2)
        +pow((x*ey)/(2*sqrt(y*y*y)),2);
    return sqrt(
        Contain
    );
}
```

**Mathematical Formula**:
$$\sigma_z = \sqrt{\left(\frac{\sigma_x}{\sqrt{y}}\right)^2 + \left(\frac{x\sigma_y}{2y^{3/2}}\right)^2}$$

**Used for**: V_n PtDiff = V_nΔ / √(V_nΔ_ref)

**Location in Code**: 
- [Process_FourierFit.cxx](Process_FourierFit.cxx#L394): Line ~398
- [Process_TemplateFit.cxx](Process_TemplateFit.cxx#L479): Line ~486

---

#### 3. Error_vNL - Flow Harmonic Error  

```cpp
// Line 28-36
double Error_vNL(double N, double eN, double D1, double eD1){
    //vNL = Numerator/sqrt(Denominator1)
    double err = sqrt(
        pow(eN/sqrt(D1),2)
        +pow(N*eD1/(2*pow(D1,3./2.)),2)
    );
    if(err<=0||err>5)printf("Warning: err for vNL is %f\n",err);
    return err;
}
```

**Mathematical Formula**:
$$\sigma_v = \sqrt{\left(\frac{\sigma_N}{\sqrt{D}}\right)^2 + \left(\frac{N\sigma_D}{2D^{3/2}}\right)^2}$$

---

## Formulas Summary

### Chi² Calculation

**Standard Formula**:
$$\chi^2 = \sum_{i=1}^{N_{bins}} \frac{(y_i^{data} - y_i^{fit})^2}{\sigma_i^2}$$

**Reduced Chi²**:
$$\frac{\chi^2}{ndf} = \frac{\chi^2}{N_{bins} - N_{params}}$$

**For all methods**: N_params = 5

---

### Fourier Fit Function

$$f(\Delta\varphi) = a_0 + 2\sum_{n=1}^{4} a_n \cos(n\Delta\varphi)$$

**Parameters**: [a₀, a₁, a₂, a₃, a₄]

**Flow extraction**:
$$v_n^2 = \frac{a_n}{a_0} \quad \Rightarrow \quad v_n = \sqrt{\frac{a_n}{a_0}}$$

---

### Template Fit Function

$$C_{HM}(\Delta\varphi) = F \cdot C_{LM}(\Delta\varphi) + G \cdot \left[1 + 2\sum_{n=2,3,4} v_{n,1}^2 \cos(n\Delta\varphi)\right]$$

where:
- C_HM: High multiplicity correlation (data)
- C_LM: Low multiplicity correlation (template)
- F: Template normalization factor
- G: Flow component amplitude
- v_{n,1}²: Squared flow harmonics

**Parameters**: [F, G, v²₁, v²₃, v²₄]

---

### Systematic Uncertainty Formula

**Combined Error**:
$$\sigma_{total} = \text{LR\_ERROR\_SCALE} \times \sqrt{\sigma_{stat}^2 + (\text{LR\_SYS\_FRAC} \times |y_i|)^2}$$

**With Floor**:
$$\sigma_{total} = \max\left(\sigma_{above}, \text{LR\_ERROR\_FLOOR}\right)$$

**Quadrature Addition**:
$$\sigma_{combined}^2 = \sigma_{statistical}^2 + \sigma_{systematic}^2$$

---

## Environment Variables Reference

### Template Fit Systematics Control

| Variable | Type | Default | Range | Purpose |
|----------|------|---------|-------|---------|
| **LR_ERROR_SCALE** | double | 1.0 | > 0 | Multiplicative factor on all errors |
| **LR_SYS_FRAC** | double | 0.0 | ≥ 0 | Fractional systematic (0.01 = 1%) |
| **LR_ERROR_FLOOR** | double | 1e-6 | ≥ 0 | Minimum error per bin |

### Usage Examples

```bash
# Example 1: Default (no systematics)
root -l -b -q Process_TemplateFit.cxx

# Example 2: Add 1% systematic
export LR_SYS_FRAC=0.01
root -l -b -q Process_TemplateFit.cxx

# Example 3: Add 2% systematic + 10% error scaling
export LR_SYS_FRAC=0.02
export LR_ERROR_SCALE=1.1
root -l -b -q Process_TemplateFit.cxx

# Example 4: Conservative analysis
export LR_SYS_FRAC=0.05
export LR_ERROR_SCALE=1.2
export LR_ERROR_FLOOR=0.001
root -l -b -q Process_TemplateFit.cxx

# Reset to defaults
unset LR_ERROR_SCALE LR_SYS_FRAC LR_ERROR_FLOOR
```

### Output Messages

When systematics are enabled, you'll see:
```
[FitDiag] Applying runtime error regularization: scale=1.1 sysFrac=0.02 floor=1e-06
```

When disabled (default):
```
(No message - systematics not applied)
```

---

## Thresholds and Quality Criteria

### Chi²/ndf Rejection Thresholds

| Method | Threshold | Action if Exceeded |
|--------|-----------|-------------------|
| **Fourier Fit** | χ²/ndf > 2000 | Set all params to -1, errors to 10, return |
| **Template Fit** | χ²/ndf > 30000 | Set all params to -1, errors to 10, return |

### Interpretation Guidelines

| χ²/ndf | Quality | Interpretation |
|--------|---------|----------------|
| ~1.0 | Excellent | Data perfectly described by model |
| 0.5-2.0 | Good | Normal statistical fluctuations |
| 2-5 | Acceptable | Some model-data tension |
| 5-10 | Marginal | Significant discrepancies |
| 10-100 | Poor | Major issues (but not rejected) |
| >2000 (Fourier) | Failed | REJECTED |
| >30000 (Template) | Failed | REJECTED |

**Note**: Template fit has higher threshold because of more complex model

---

## Data Flow: Histogram → V_n

```
Step 1: Histogram Bin Errors
    hm->GetBinError(i)  →  σ_stat = √N
                ↓
Step 2: Add Systematics (optional)
    if (LR_SYS_FRAC > 0):
        σ_sys = LR_SYS_FRAC × |bin_content|
        σ_total = √(σ_stat² + σ_sys²)
    else:
        σ_total = σ_stat
                ↓
Step 3: Chi² Fit
    χ² = Σ[(data - fit)²/σ_total²]
    Fit parameters: [a₀±σ_a₀, a₁±σ_a₁, ..., a₄±σ_a₄]
                ↓
Step 4: Extract v_n²
    v_n² = a_n / a_0
    σ_v_n² = Error_Ratio(a_n, σ_a_n, a_0, σ_a_0, rho)
                ↓
Step 5: Convert to v_n
    v_n = √(v_n²)
    σ_v_n = ½ × σ_v_n² / v_n
                ↓
Step 6: PtDiff (if applicable)
    V_n = V_nΔ / √(V_nΔ_ref)
    σ_V_n = Error_Ratio_sqrtY(V_nΔ, σ_V_nΔ, V_nΔ_ref, σ_V_nΔ_ref)
                ↓
Final Result: V_n ± σ_total
    (includes systematics if added in Step 2)
```

---

## Output Files and Plots

### Generated Files

**From Fourier Fit**:
```
./FourierFit/
├── Vn_[dataset]_[mult_range].root       # Output ROOT file with V_n
└── PDFs/
    └── FourierFit_[dataset]_[mult_range].pdf  # Fit quality plot
```

**From Template Fit**:
```
./TemplateFit/
├── Vn_Template_[dataset]_[mult_range].root
└── PDFs/
    └── TemplateFit_[dataset]_[mult_range].pdf
```

### What's in the Plots

Both methods generate plots showing:
1. **Top panel**: Data points with fit curve and χ²/ndf value
2. **Bottom panel**: Residuals or pull distribution
3. **Text labels**: 
   - χ²/ndf = X.X/N = Y.YY
   - Extracted V_n values
   - Fit parameters

---

## Common Issues and Debugging

### Problem: χ²/ndf is too large (but below threshold)

**Possible causes**:
- Template from wrong centrality bin
- Insufficient statistics
- Wrong fit range
- Model doesn't match data

**Solutions**:
- Check template selection
- Try different fit ranges
- Add systematic uncertainties (increase LR_SYS_FRAC)
- Inspect residual plots

---

### Problem: Fit fails (exceeds threshold)

**Diagnosis**:
```cpp
if (chi2ndf > threshold) {
    // All parameters set to -1
    // Check output ROOT file: v2 = -1 means failed fit
}
```

**Solutions**:
- Verify input data quality
- Check for empty bins
- Try different initial parameters
- Consider alternative fit method

---

### Problem: Systematic uncertainties not applied

**Symptoms**:
- No "[FitDiag]" message in output
- Errors unchanged from √N

**Check**:
```bash
# Verify environment variables are set
echo $LR_SYS_FRAC
echo $LR_ERROR_SCALE
echo $LR_ERROR_FLOOR

# If empty, systematics are disabled
```

**Solution**:
```bash
# Set before running
export LR_SYS_FRAC=0.01
root -l -b -q Process_TemplateFit.cxx
```

---

### Problem: Errors seem too small

**Check**:
1. Bin content vs bin error ratio
2. Bootstrap errors (if using 3times2PC)
3. Whether systematics are included

**Compare**:
```
Without syst: V_n = 0.045 ± 0.001
With 2% syst: V_n = 0.045 ± 0.002
```

---

## Quick Command Reference

```bash
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Full analysis workflow
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# Step 1: Extract correlations
root -l -b -q Process_dPhidEta.cxx

# Step 2: Create bootstrap samples
root -l -b -q Process_CreateBootstrapSample.cxx

# Step 3a: Fourier fit (no systematics)
root -l -b -q Process_FourierFit.cxx

# Step 3b: Template fit (default: no systematics)
root -l -b -q Process_TemplateFit.cxx

# Step 3c: Template fit WITH 1.5% systematics
export LR_SYS_FRAC=0.015
root -l -b -q Process_TemplateFit.cxx
unset LR_SYS_FRAC

# Optional: 3times2PC method
root -l -b -q Process_3times2PC.cxx

#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Run with different systematic settings
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# Conservative: 3% systematic
export LR_SYS_FRAC=0.03
root -l -b -q Process_TemplateFit.cxx

# Very conservative: 5% systematic + 20% error scale
export LR_SYS_FRAC=0.05
export LR_ERROR_SCALE=1.2
root -l -b -q Process_TemplateFit.cxx

# Reset
unset LR_SYS_FRAC LR_ERROR_SCALE

#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Check output
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# View fit quality plots
evince ./TemplateFit/PDFs/TemplateFit_*.pdf &

# Check ROOT file
root -l ./TemplateFit/Vn_Template_*.root
# Inside ROOT:
# .ls
# hV2->Draw()
# hV2->Print("all")
```

---

## Contact Information

**Framework Authors**:
- Main framework: Zhiyong Lu (zhiyong.lu@cern.ch)
- Template fitter: Vytautas Vislavicius (vytautas.vislavicius@cern.ch)

**For questions about**:
- Chi² calculation: See this guide
- Systematic uncertainties: Check environment variables section
- Error propagation: See [include/ErrorPropagation.h](include/ErrorPropagation.h)
- Fit quality: Check χ²/ndf thresholds section

---

## Related Documentation

1. **CHI2_CALCULATION_OVERVIEW.md** - Complete technical documentation
2. **CHI2_VISUAL_GUIDE.md** - Step-by-step visual workflows
3. **README.md** - General usage instructions

