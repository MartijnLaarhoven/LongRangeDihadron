# Chi² Calculation Overview - Long Range Dihadron Analysis

## Table of Contents
1. [Executive Summary](#executive-summary)
2. [Chi² Calculation Methods](#chi-calculation-methods)
3. [Step-by-Step Chi² Calculation](#step-by-step-chi-calculation)
4. [Systematic Uncertainties](#systematic-uncertainties)
5. [Error Propagation](#error-propagation)
6. [Summary Tables](#summary-tables)

---

## Executive Summary

### What is Chi²?
Chi-squared (χ²) is a goodness-of-fit metric that quantifies how well a fitted function matches the data:

$$\chi^2 = \sum_{i=1}^{N_{bins}} \frac{(y_i^{data} - y_i^{fit})^2}{\sigma_i^2}$$

where:
- $y_i^{data}$ = data bin content
- $y_i^{fit}$ = fitted function value at bin center
- $\sigma_i$ = uncertainty on bin content
- $N_{bins}$ = number of bins

**Reduced Chi²**: $\chi^2/ndf$ where $ndf = N_{bins} - N_{params}$

A good fit typically has $\chi^2/ndf \approx 1$.

### Where Chi² is Calculated

| Module | Function | Line | Has Systematics? |
|--------|----------|------|------------------|
| **Fourier Fit** | `RooFourierFit()` | ~681 | ❌ No |
| **Fourier Fit** | `PlotFitting()` | ~838 | ❌ No |
| **Template Fit** | `RooTemplateFit()` | ~500+ | ✅ Yes (optional) |
| **Template Fit** | `PlotFitting()` | ~1088 | ✅ Yes (optional) |
| **Generic Fitter** | `TemplateFitter::Fit()` | ~89 | ✅ Yes (optional) |

---

## Chi² Calculation Methods

### Method 1: Fourier Fit (Process_FourierFit.cxx)

Two implementations exist in this file:

#### 1A. RooFit Native Chi² (Line ~681)

**Location**: `void RooFourierFit(TH1 *hist, ...)`

```cpp
// After fit completes
double chi2 = fitFunc->GetChisquare();
double ndf = fitFunc->GetNDF();
double chi2_ndf = (ndf > 0) ? chi2 / ndf : 0;
std::cout << "Chi2/NDF = " << chi2_ndf << std::endl;

// Quality check
if (chi2_ndf > 2000) {
    std::cout << "WARNING: Chi2/NDF > 20, fit may be inaccurate!" << std::endl;
    // Set all parameters to -1 (failed fit indicator)
    return;
}
```

**Fit Function**:
$$f(x) = a_0 + 2a_1\cos(x) + 2a_2\cos(2x) + 2a_3\cos(3x) + 2a_4\cos(4x)$$

**Steps**:
1. ROOT's `TF1::Fit()` performs least-squares minimization
2. `GetChisquare()` returns the χ² value from the fit
3. `GetNDF()` returns degrees of freedom (nBins - 5)
4. Reject if χ²/ndf > 2000 (extremely poor fit)

**Error Model**: Bin errors are taken as-is from histogram (typically $\sqrt{N}$ from Poisson statistics)

---

#### 1B. Manual Chi² Calculation (Line ~838)

**Location**: `void PlotFitting(TH1 *hm, ...)`

```cpp
// =============== Calculate chi2/ndf ===============
double chi2 = 0.0;
int nBins = hm->GetNbinsX();
int nParams = 5; // a0, a1, a2, a3, a4
int ndf = nBins - nParams;

for (int i = 1; i <= nBins; i++) {
    double data = hm->GetBinContent(i);
    double error = hm->GetBinError(i);
    double x = hm->GetBinCenter(i);
    double fit = fit_p1m->Eval(x);
    
    if (error > 0) { // Skip bins with zero error
        double residual = data - fit;
        chi2 += (residual * residual) / (error * error);
    }
}
double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;
```

**Step-by-Step**:
1. Loop over all bins in the histogram
2. For each bin:
   - Get data value: `data = hm->GetBinContent(i)`
   - Get uncertainty: `error = hm->GetBinError(i)`
   - Get bin center position: `x = hm->GetBinCenter(i)`
   - Evaluate fitted function: `fit = fit_p1m->Eval(x)`
   - Calculate contribution: $\chi^2_i = \frac{(data - fit)^2}{error^2}$
3. Sum all contributions: $\chi^2 = \sum_i \chi^2_i$
4. Calculate degrees of freedom: $ndf = N_{bins} - 5$
5. Calculate reduced chi²: $\chi^2/ndf$

---

### Method 2: Template Fit (Process_TemplateFit.cxx)

#### 2A. RooFit Template Chi² (Line ~500+)

**Location**: `void RooTemplateFit(TH1 *lm, TH1 *hm, ...)`

This is more complex because it allows **systematic uncertainties** to be added before fitting.

**IMPORTANT**: Systematic uncertainties are controlled by environment variables:
- `LR_ERROR_SCALE` - multiplicative factor on all errors (default: 1.0)
- `LR_SYS_FRAC` - fractional systematic uncertainty (default: 0.0)
- `LR_ERROR_FLOOR` - minimum error per bin (default: 1e-6)

**Template Fit Model**:
$$C_{high\,mult}(\Delta\varphi) = F \cdot C_{low\,mult}(\Delta\varphi) + G \cdot [1 + 2v_{2,1}^2\cos(2\Delta\varphi) + 2v_{3,1}^2\cos(3\Delta\varphi) + 2v_{4,1}^2\cos(4\Delta\varphi)]$$

where:
- $C_{high\,mult}$: High multiplicity correlation
- $C_{low\,mult}$: Low multiplicity correlation (template)
- $F$: Template normalization factor
- $G$: Flow component amplitude
- $v_{n,1}^2$: Squared flow harmonics

---

#### 2B. Manual Chi² for Template Fit (Line ~1088)

**Location**: `void PlotFitting(TH1 *lm, TH1 *hm, ...)`

```cpp
TH1D* hsubtract = (TH1D*)hm->Clone(Form("subtract"));
hsubtract->Add(lm, -F);  // Subtract template: data - F*template

// =============== Calculate chi2/ndf ===============
double chi2 = 0.0;
int nBins = hsubtract->GetNbinsX();
int nParams = 5; // F, G, v21, v31, v41
int ndf = nBins - nParams;

for (int i = 1; i <= nBins; i++) {
    double data = hsubtract->GetBinContent(i);
    double error = hsubtract->GetBinError(i);
    double x = hsubtract->GetBinCenter(i);
    double fit = fit_p1m->Eval(x);
    
    if (error > 0) {
        double residual = data - fit;
        chi2 += (residual * residual) / (error * error);
    }
}
double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;
```

**Key Difference**: This fits the **background-subtracted** histogram, not the raw data.

---

## Step-by-Step Chi² Calculation

### Complete Workflow for Template Fit (with Systematics)

#### STEP 1: Load Input Histograms

```cpp
TH1* lm = GetHistogram("low_multiplicity");   // Template
TH1* hm = GetHistogram("high_multiplicity");  // Data
```

---

#### STEP 2: Apply Systematic Uncertainties (OPTIONAL)

**Code Location**: Process_TemplateFit.cxx, lines 748-776

```cpp
// Read environment variables
const char* es = getenv("LR_ERROR_SCALE");
const char* sf = getenv("LR_SYS_FRAC");
const char* ef = getenv("LR_ERROR_FLOOR");

double envScale = es ? atof(es) : 1.0;
double sysFrac = sf ? atof(sf) : 0.0;
double errFloor = ef ? atof(ef) : 1e-6;
```

**For each bin**:

```cpp
for (int ib = 1; ib <= nb; ++ib) {
    // Get original statistical error
    double err = hm->GetBinError(ib);      // σ_stat = √N
    double val = hm->GetBinContent(ib);    // bin content
    
    // If error is invalid, estimate from content
    if (!TMath::Finite(err) || err <= 0) {
        err = (val > 0) ? sqrt(val) : 1.0;
    }
    
    // Add systematic in quadrature
    if (sysFrac > 0.0) {
        double sys = fabs(val) * sysFrac;              // σ_sys = |content| × sysFrac
        err = sqrt(err*err + sys*sys);                 // σ_total = √(σ_stat² + σ_sys²)
    }
    
    // Apply scaling factor
    err *= envScale;
    
    // Apply error floor
    if (!TMath::Finite(err) || err < errFloor) {
        err = errFloor;
    }
    
    // Update bin error
    hm->SetBinError(ib, err);
}
```

**Mathematical Formula**:
$$\sigma_{total}(i) = \max\left( \text{LR\_ERROR\_SCALE} \times \sqrt{\sigma_{stat}^2(i) + \left(\text{LR\_SYS\_FRAC} \times |y_i|\right)^2}, \, \text{LR\_ERROR\_FLOOR} \right)$$

**Example with Numbers**:

| Bin | Content | $\sigma_{stat}$ | LR_SYS_FRAC | $\sigma_{sys}$ | $\sigma_{total}$ |
|-----|---------|-----------------|-------------|----------------|------------------|
| 1 | 100 | 10 | 0.01 (1%) | 1.0 | √(100+1) = 10.05 |
| 2 | 1000 | 31.6 | 0.01 (1%) | 10.0 | √(1000+100) = 33.2 |
| 3 | 10 | 3.16 | 0.01 (1%) | 0.1 | √(10+0.01) = 3.16 |

**Impact**: Larger bins get relatively more systematic uncertainty, reducing their weight in the fit.

---

#### STEP 3: Normalize Template to Data

```cpp
// Scale template to match data mean
double scale_lm_to_hm = hm_mean / lm_mean;
lm->Scale(scale_lm_to_hm);
```

**Why**: Makes initial parameter $F \approx 1$, improving fit stability.

---

#### STEP 4: Estimate Initial Parameters

```cpp
// Weighted linear regression: hm ≈ Fa * lm + Ga
double S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
for (int ib = 1; ib <= nbins; ++ib) {
    double x = lm->GetBinContent(ib);
    double y = hm->GetBinContent(ib);
    double w = 1.0 / (error * error);  // Weight = 1/σ²
    
    S += w;
    Sx += w * x;
    Sy += w * y;
    Sxx += w * x * x;
    Sxy += w * x * y;
}

// Solve: Fa = (S*Sxy - Sx*Sy) / (S*Sxx - Sx*Sx)
//        Ga = (Sy - Fa*Sx) / S
double Fa = (S * Sxy - Sx * Sy) / (S * Sxx - Sx * Sx);
double Ga = (Sy - Fa * Sx) / S;
```

---

#### STEP 5: Perform RooFit Chi² Fit

```cpp
// Create RooFit objects
RooRealVar x("x", "x", xmin, xmax);
RooDataHist dsig("dsig", "dsig", x, Import(*hm));

// Define fit function with RooFit
RooRealVar F("F", "F", Fa_init, 0.1*Fa_init, 10*Fa_init);
RooRealVar G("G", "G", Ga_init, 0.1*Ga_init, 10*Ga_init);
RooRealVar v21("v21", "v21", 0.01, 0, 0.5);
RooRealVar v31("v31", "v31", 0.01, 0, 0.5);
RooRealVar v41("v41", "v41", 0.01, 0, 0.5);

// Fit function
RooAbsReal* fitFunc = /* define template + flow function */

// Perform chi² fit
fitFunc->chi2FitTo(dsig, 
                   RooFit::PrintLevel(-1), 
                   RooFit::SumW2Error(kTRUE),  // Use SumW2 weights
                   RooFit::Warnings(kFALSE));
```

---

#### STEP 6: Calculate Chi²/ndf

```cpp
RooAbsReal* chi2_o_ndf = fitFunc->createChi2(dsig, 
                                              Range("cutrange"), 
                                              Extended(true), 
                                              DataError(RooAbsData::SumW2));

double chi2 = chi2_o_ndf->getVal();
int ndf = dataH->GetNbinsX() - 5;  // nBins - nParams
double chi2ndf = chi2 / ndf;

std::cout << "chi2: " << chi2 << std::endl;
std::cout << "ndf: " << ndf << ", chi2/ndf: " << chi2ndf << std::endl;
```

**RooFit Chi² Formula** (SumW2 mode):
$$\chi^2 = \sum_{i=1}^{N_{bins}} \frac{(y_i - f_i)^2}{\sigma_i^2}$$

where $\sigma_i$ is the **modified** bin error (includes systematics if LR_SYS_FRAC > 0).

---

#### STEP 7: Quality Check

```cpp
if (chi2ndf > 30000) {  // Template fit threshold
    // Reject fit - set all parameters to -1
    for(Int_t i=0; i<fParList->GetEntries(); i++) {
        ((RooRealVar*)fParList->At(i))->setVal(-1);
        ((RooRealVar*)fParList->At(i))->setError(10);
    }
    return 1;
}
```

**Thresholds**:
- Fourier Fit: χ²/ndf > 2000 → REJECT
- Template Fit: χ²/ndf > 30000 → REJECT

---

#### STEP 8: Extract Fit Parameters

```cpp
// Get fitted parameters
double F_fitted = F.getVal();
double F_error = F.getError();
double G_fitted = G.getVal();
double G_error = G.getError();
double v21_fitted = v21.getVal();
double v21_error = v21.getError();
// ... etc for v31, v41

// Store in output structure
vnResult->v2 = sqrt(v21_fitted);  // Convert v²₂ to v₂
vnResult->v2_err = 0.5 * v21_error / sqrt(v21_fitted);  // Error propagation
```

**Note**: The parameter errors (`F_error`, `v21_error`, etc.) **include the effect of systematics** if they were added in Step 2.

---

## Systematic Uncertainties

### Implementation (Template Fit Only)

Systematic uncertainties are **NOT hardcoded** in the analysis. They are controlled at runtime via environment variables.

#### Environment Variables

| Variable | Meaning | Default | Units |
|----------|---------|---------|-------|
| `LR_ERROR_SCALE` | Multiplicative scaling of all errors | 1.0 | dimensionless |
| `LR_SYS_FRAC` | Fractional systematic uncertainty | 0.0 | fraction (0.01 = 1%) |
| `LR_ERROR_FLOOR` | Minimum allowed error per bin | 1e-6 | absolute |

#### How to Use

```bash
# Default: No systematics
root -l -b -q Process_TemplateFit.cxx

# Add 1% systematic uncertainty
export LR_SYS_FRAC=0.01
root -l -b -q Process_TemplateFit.cxx

# Add 2% systematic + 10% error scaling
export LR_SYS_FRAC=0.02
export LR_ERROR_SCALE=1.1
root -l -b -q Process_TemplateFit.cxx

# Reset to defaults
unset LR_ERROR_SCALE LR_SYS_FRAC LR_ERROR_FLOOR
```

#### Output Indication

When systematics are applied, you'll see:
```
[FitDiag] Applying runtime error regularization: scale=1.0 sysFrac=0.01 floor=1e-06
```

---

### What Is NOT Implemented

The following systematic uncertainties are **NOT** handled by the current framework:

❌ **Detector Systematics**:
- Track reconstruction efficiency variations
- PID selection efficiency
- Detector acceptance effects
- Dead channel corrections

❌ **Method Systematics**:
- Variation of template shape
- Fit range dependence
- Background subtraction method variations
- Alternative functional forms

❌ **Correlated Systematics**:
- Covariance matrices between bins
- Common-mode uncertainties
- Systematic eigenvector decomposition

❌ **Physics Systematics**:
- MC generator dependence
- Resonance contribution variations
- Flow factorization breaking

**These must be studied separately** and combined with the statistical uncertainties manually.

---

## Error Propagation

### Error Propagation Functions

**Location**: `include/ErrorPropagation.h`

The code includes sophisticated error propagation for derived quantities:

#### 1. Simple Ratio: z = x/y

```cpp
double Error_Ratio(double x, double ex, double y, double ey, double rho)
{
    // rho = correlation coefficient
    // Cov(x,y) = rho * ex * ey
    
    double Contain = pow(ex/y, 2)
                   + pow(x*ey/(y*y), 2)
                   + 2*(1./y)*(-1.*x/(y*y))*rho*ex*ey;
    
    return sqrt(Contain);
}
```

**Formula**:
$$\sigma_z^2 = \left(\frac{\sigma_x}{y}\right)^2 + \left(\frac{x\sigma_y}{y^2}\right)^2 + 2\frac{\partial z}{\partial x}\frac{\partial z}{\partial y}\rho\sigma_x\sigma_y$$

---

#### 2. Ratio with Square Root: z = x/√y

```cpp
double Error_Ratio_sqrtY(double x, double ex, double y, double ey)
{
    // Used for: V_n = V_nΔ / √(V_nΔ_ref)
    
    return sqrt(
        pow(ex/sqrt(y), 2)
        + pow((x*ey)/(2*sqrt(y*y*y)), 2)
    );
}
```

**Formula**:
$$\sigma_z = \sqrt{\left(\frac{\sigma_x}{\sqrt{y}}\right)^2 + \left(\frac{x\sigma_y}{2y^{3/2}}\right)^2}$$

**Used for**: Converting $V_{n\Delta}$ to $V_n$ in pT-differential measurements.

---

#### 3. Three-Particle Correlation Error

```cpp
double Error_SCklm(double N_klm, double eN_klm, 
                   double N_kl, double eN_kl, 
                   double N_km, double eN_km, 
                   double N_lm, double eN_lm,
                   double N_k, double eN_k, 
                   double N_l, double eN_l, 
                   double N_m, double eN_m)
{
    // SC(k,l,m) = N_klm - N_kl*N_m - N_km*N_l - N_lm*N_k + 2*N_k*N_l*N_m
    
    double err = sqrt(
        pow(eN_klm, 2) + pow(eN_kl*N_m, 2) + pow(eN_km*N_l, 2) + pow(eN_lm*N_k, 2)
        + pow((2*N_l*N_m - N_lm)*eN_k, 2)
        + pow((2*N_k*N_m - N_km)*eN_l, 2)
        + pow((2*N_k*N_l - N_kl)*eN_m, 2)
    );
    
    return err;
}
```

**Used for**: 3times2PC method error propagation.

---

### Error Flow: From Histogram to V_n

```
Histogram Bin Errors (√N)
         ↓
    [Optional: Add Systematics]
    σ_total = √(σ_stat² + σ_sys²)
         ↓
    Chi² Fit → Fit Parameter Errors
    (F±σ_F, G±σ_G, v²₁±σ_v²₁, ...)
         ↓
    Error Propagation
    v_n = √(v²_n)
    σ_vn = ½ σ_v²n / √(v²_n)
         ↓
    [If PtDiff: Divide by Reference]
    V_n = V_nΔ / √(V_nΔ_ref)
    σ_Vn = Error_Ratio_sqrtY(...)
         ↓
    Final V_n ± Statistical Error
    (includes systematic if added in step 1)
```

---

## Summary Tables

### Chi² Implementations

| File | Function | Line | Method | Has Syst? | Threshold |
|------|----------|------|--------|-----------|-----------|
| Process_FourierFit.cxx | RooFourierFit() | ~681 | RooFit native | ❌ | χ²/ndf > 2000 |
| Process_FourierFit.cxx | PlotFitting() | ~838 | Manual loop | ❌ | (display only) |
| Process_TemplateFit.cxx | RooTemplateFit() | ~500+ | RooFit native | ✅* | χ²/ndf > 30000 |
| Process_TemplateFit.cxx | PlotFitting() | ~1088 | Manual loop | ✅* | (display only) |
| include/TemplateFitter.cxx | Fit() | ~89 | RooFit native | ✅* | χ²/ndf > 30000 |

*Systematics via environment variables (disabled by default)

---

### Degree of Freedom Calculation

All methods use the same formula:

$$ndf = N_{bins} - N_{params}$$

| Method | N_params | Parameters |
|--------|----------|------------|
| Fourier Fit | 5 | a₀, a₁, a₂, a₃, a₄ |
| Template Fit | 5 | F, G, v²₁, v²₃, v²₄ |

---

### Error Models

| Module | Error Model | Formula |
|--------|-------------|---------|
| All modules | Poisson (default) | σᵢ = √Nᵢ |
| Template Fit | + Systematic (optional) | σᵢ = √(Nᵢ + (frac × Nᵢ)²) |
| RooFit | SumW2 | σᵢ² = Σwᵢ² |

---

## Recommendations

### For Your Analysis

1. **Document systematics**: If using LR_SYS_FRAC, document the value and justification
2. **Study fit quality**: Monitor χ²/ndf distributions to identify problematic fits
3. **Cross-check methods**: Compare Fourier Fit vs Template Fit results
4. **Extend systematics**: Consider implementing similar framework for Fourier Fit

### For Publication

1. **Separate statistical and systematic errors**: 
   - Statistical: from fit (without LR_SYS_FRAC)
   - Systematic: from variations (method, cuts, etc.)
   
2. **Quote both**: V_n = X.XXX ± 0.YYY (stat) ± 0.ZZZ (syst)

3. **Document all systematic sources**:
   - Template shape variations
   - Fit range variations
   - Selection cut variations
   - etc.

4. **Provide covariance matrices** for multi-dimensional fits if available

---

## Quick Reference: Enabling Systematics

```bash
# No systematics (default)
root -l -b -q Process_TemplateFit.cxx

# With 1.5% systematic
export LR_SYS_FRAC=0.015
root -l -b -q Process_TemplateFit.cxx

# With 2% systematic + 5% error scale
export LR_SYS_FRAC=0.02
export LR_ERROR_SCALE=1.05
root -l -b -q Process_TemplateFit.cxx

# Reset
unset LR_SYS_FRAC LR_ERROR_SCALE LR_ERROR_FLOOR
```

---

## Contact

For questions about this framework:
- Original author: Zhiyong Lu (zhiyong.lu@cern.ch)
- Template Fitter: Vytautas Vislavicius (vytautas.vislavicius@cern.ch)

