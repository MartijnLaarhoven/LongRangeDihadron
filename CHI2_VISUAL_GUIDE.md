# Chi² Calculation - Visual Step-by-Step Guide

## Table of Contents
1. [Overall Analysis Flow](#overall-analysis-flow)
2. [Step-by-Step: Fourier Fit Chi²](#step-by-step-fourier-fit-chi)
3. [Step-by-Step: Template Fit Chi²](#step-by-step-template-fit-chi)
4. [Systematic Uncertainty Addition](#systematic-uncertainty-addition)
5. [Numerical Examples](#numerical-examples)

---

## Overall Analysis Flow

```
┌─────────────────────────────────────────────────────────────────┐
│                     INPUT: ROOT Files                           │
│                  (AnalysisResults.root)                         │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│               STEP 1: Process_dPhidEta.cxx                      │
│                                                                  │
│  • Extract THnSparse histograms                                │
│  • Project to 2D: (Δφ, Δη)                                    │
│  • Separate by multiplicity/centrality                         │
│  • Create same-event and mixed-event correlations             │
│                                                                  │
│  OUTPUT: Correlation2D histograms                              │
│          (dPhidEtaSE, dPhidEtaME)                              │
└──────────────────────────┬──────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│      STEP 2: Process_CreateBootstrapSample.cxx                 │
│                                                                  │
│  • Create N bootstrap samples (typically 100)                  │
│  • Resample with replacement                                   │
│  • Used for robust error estimation                            │
│                                                                  │
│  OUTPUT: Bootstrap samples in separate folders                │
└──────────────────────────┬──────────────────────────────────────┘
                           │
           ┌───────────────┴───────────────┐
           │                               │
           ▼                               ▼
┌────────────────────────┐      ┌────────────────────────┐
│  STEP 3A: Fourier Fit  │      │  STEP 3B: Template Fit │
│                        │      │                        │
│  Process_FourierFit    │      │  Process_TemplateFit   │
│  .cxx                  │      │  .cxx                  │
│                        │      │                        │
│  ❌ NO Systematics     │      │  ✅ Systematics        │
│    (statistical only)  │      │    (optional)          │
└────────┬───────────────┘      └────────┬───────────────┘
         │                               │
         │      ╔════════════════════════╤═══════════╗
         │      ║                 Chi² FIT            ║
         │      ║                                     ║
         └──────╫──> Extract V_n harmonics           ║
                ║     (v₂, v₃, v₄)                   ║
                ║                                     ║
                ║     Calculate uncertainties         ║
                ║     (statistical ± systematic)      ║
                ╚═════════════════════════════════════╝
                           │
                           ▼
                ┌──────────────────────┐
                │   FINAL RESULTS      │
                │                      │
                │   V_n vs pT          │
                │   V_n vs multiplicity│
                │   With errors        │
                │   Chi²/ndf plots     │
                └──────────────────────┘
```

---

## Step-by-Step: Fourier Fit Chi²

### Visual Workflow

```
┌─────────────────────────────────────────────────────────────┐
│ STEP 1: Load Histogram                                      │
│                                                              │
│   TH1* hm = GetHistogram("correlation_vs_dphi");           │
│                                                              │
│   Data points: (φ₁, y₁±σ₁), (φ₂, y₂±σ₂), ..., (φₙ, yₙ±σₙ) │
│                                                              │
│   where σᵢ = hm->GetBinError(i)  [typically √N]            │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│ STEP 2: Define Fit Function                                 │
│                                                              │
│   f(φ) = a₀ + 2a₁cos(φ) + 2a₂cos(2φ) + 2a₃cos(3φ) + 2a₄cos(4φ) │
│                                                              │
│   Parameters to fit: [a₀, a₁, a₂, a₃, a₄]                  │
│                                                              │
│   Initial guesses:                                          │
│     a₀ = histogram mean                                     │
│     a₁ = 0.001, a₂ = 0.001, a₃ = 0.001, a₄ = 0.001       │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│ STEP 3: Perform Least-Squares Fit                          │
│                                                              │
│   Minimize:                                                  │
│                                                              │
│   χ² = Σᵢ [ (yᵢ - f(φᵢ))² / σᵢ² ]                         │
│                                                              │
│   ROOT's TF1::Fit() does this automatically using:         │
│   - MINUIT minimizer                                        │
│   - Iterative optimization                                  │
│   - Numerical derivatives                                   │
│                                                              │
│   Result: Best-fit parameters [a₀, a₁, a₂, a₃, a₄]        │
│           with uncertainties [σ_a₀, σ_a₁, ..., σ_a₄]      │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│ STEP 4: Calculate Chi² from Fit                            │
│                                                              │
│   Method A (RooFit native):                                 │
│   ────────────────────────                                  │
│   chi2 = fitFunc->GetChisquare();                          │
│   ndf = fitFunc->GetNDF();                                 │
│                                                              │
│   Method B (Manual calculation):                            │
│   ────────────────────────────                              │
│   chi2 = 0.0;                                              │
│   for (i = 1; i <= nBins; i++) {                          │
│       data = hm->GetBinContent(i);                         │
│       fit = fitFunc->Eval(hm->GetBinCenter(i));           │
│       error = hm->GetBinError(i);                          │
│       if (error > 0) {                                     │
│           chi2 += pow((data - fit) / error, 2);           │
│       }                                                     │
│   }                                                         │
│   ndf = nBins - 5;                                         │
│   chi2ndf = chi2 / ndf;                                    │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│ STEP 5: Quality Check                                       │
│                                                              │
│   if (chi2/ndf > 2000) {                                   │
│       // REJECT FIT - Too poor quality                     │
│       Set all parameters to -1                             │
│       Set all errors to 10                                 │
│       return;                                              │
│   }                                                         │
│                                                              │
│   Good fits typically have: 0.5 < χ²/ndf < 2.0            │
│   Acceptable: χ²/ndf < 10                                  │
│   Poor but kept: 10 < χ²/ndf < 2000                       │
│   Rejected: χ²/ndf > 2000                                  │
└─────────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────────┐
│ STEP 6: Extract Flow Coefficients                          │
│                                                              │
│   From fitted parameters, calculate:                        │
│                                                              │
│   v₁² = a₁ / a₀    →    v₁ = √(a₁/a₀)                     │
│   v₂² = a₂ / a₀    →    v₂ = √(a₂/a₀)                     │
│   v₃² = a₃ / a₀    →    v₃ = √(a₃/a₀)                     │
│   v₄² = a₄ / a₀    →    v₄ = √(a₄/a₀)                     │
│                                                              │
│   Errors via Error_Ratio():                                 │
│   σ(v_n²) = Error_Ratio(aₙ, σ_aₙ, a₀, σ_a₀, 0)           │
│   σ(v_n) = ½ × σ(v_n²) / v_n                              │
└─────────────────────────────────────────────────────────────┘
                           ↓
                    ┌─────────────┐
                    │ FINAL V_n   │
                    │ with errors │
                    └─────────────┘
```

---

## Step-by-Step: Template Fit Chi²

### Visual Workflow (with Systematic Uncertainties)

```
┌──────────────────────────────────────────────────────────────────┐
│ STEP 1: Load Two Histograms                                     │
│                                                                   │
│   TH1* lm = GetHistogram("low_multiplicity");  // Template      │
│   TH1* hm = GetHistogram("high_multiplicity"); // Data          │
│                                                                   │
│   Both have bin errors: σᵢ = √Nᵢ                               │
└──────────────────────────────────────────────────────────────────┘
                           ↓
┌──────────────────────────────────────────────────────────────────┐
│ STEP 2: Check for Systematic Uncertainty Settings               │
│                                                                   │
│   Read environment variables:                                    │
│   ┌────────────────────────────────────────────────────────┐   │
│   │ LR_ERROR_SCALE  = getenv("LR_ERROR_SCALE")  ─┐        │   │
│   │                   default: 1.0                │        │   │
│   │                                               │        │   │
│   │ LR_SYS_FRAC     = getenv("LR_SYS_FRAC")    ──┼─ KEY!  │   │
│   │                   default: 0.0                │        │   │
│   │                                               │        │   │
│   │ LR_ERROR_FLOOR  = getenv("LR_ERROR_FLOOR")  ─┘        │   │
│   │                   default: 1e-6                        │   │
│   └────────────────────────────────────────────────────────┘   │
│                                                                   │
│   If LR_SYS_FRAC = 0.0  →  NO systematics (default)             │
│   If LR_SYS_FRAC > 0.0  →  ADD systematics in quadrature        │
└──────────────────────────────────────────────────────────────────┘
                           ↓
┌──────────────────────────────────────────────────────────────────┐
│ STEP 3: Apply Systematic Uncertainties (if enabled)             │
│                                                                   │
│   FOR each bin i = 1 to nBins:                                  │
│   ┌──────────────────────────────────────────────────────┐     │
│   │  (1) Get original values:                            │     │
│   │      σ_stat = hm->GetBinError(i)                     │     │
│   │      value = hm->GetBinContent(i)                    │     │
│   │                                                       │     │
│   │  (2) Calculate systematic:                           │     │
│   │      if (LR_SYS_FRAC > 0):                           │     │
│   │          σ_sys = |value| × LR_SYS_FRAC               │     │
│   │      else:                                            │     │
│   │          σ_sys = 0                                    │     │
│   │                                                       │     │
│   │  (3) Combine in quadrature:                          │     │
│   │      σ_combined = √(σ_stat² + σ_sys²)                │     │
│   │                                                       │     │
│   │  (4) Apply error scale:                              │     │
│   │      σ_combined = σ_combined × LR_ERROR_SCALE        │     │
│   │                                                       │     │
│   │  (5) Apply error floor:                              │     │
│   │      if (σ_combined < LR_ERROR_FLOOR):               │     │
│   │          σ_combined = LR_ERROR_FLOOR                 │     │
│   │                                                       │     │
│   │  (6) Update histogram:                               │     │
│   │      hm->SetBinError(i, σ_combined)                  │     │
│   └──────────────────────────────────────────────────────┘     │
│                                                                   │
│   RESULT: Modified histogram with enlarged errors                │
└──────────────────────────────────────────────────────────────────┘
                           ↓
┌──────────────────────────────────────────────────────────────────┐
│ STEP 4: Normalize Template to Data                              │
│                                                                   │
│   Calculate means:                                               │
│   hm_mean = Average of hm over all bins                         │
│   lm_mean = Average of lm over all bins                         │
│                                                                   │
│   Scale template:                                                │
│   scale = hm_mean / lm_mean                                     │
│   lm->Scale(scale)                                              │
│                                                                   │
│   PURPOSE: Makes initial parameter F ≈ 1 for better stability   │
└──────────────────────────────────────────────────────────────────┘
                           ↓
┌──────────────────────────────────────────────────────────────────┐
│ STEP 5: Estimate Initial Parameters (Linear Regression)         │
│                                                                   │
│   Model: hm ≈ F × lm + G                                         │
│                                                                   │
│   Weighted least-squares:                                        │
│   ┌──────────────────────────────────────────────────────┐     │
│   │  S   = Σᵢ wᵢ               where wᵢ = 1/σᵢ²          │     │
│   │  Sx  = Σᵢ wᵢ × lmᵢ                                   │     │
│   │  Sy  = Σᵢ wᵢ × hmᵢ                                   │     │
│   │  Sxx = Σᵢ wᵢ × lmᵢ²                                  │     │
│   │  Sxy = Σᵢ wᵢ × lmᵢ × hmᵢ                             │     │
│   │                                                       │     │
│   │  F = (S×Sxy - Sx×Sy) / (S×Sxx - Sx²)                │     │
│   │  G = (Sy - F×Sx) / S                                 │     │
│   └──────────────────────────────────────────────────────┘     │
│                                                                   │
│   Initial guesses: F_init = F, G_init = G                       │
│                    v₂ = 0.01, v₃ = 0.01, v₄ = 0.01             │
└──────────────────────────────────────────────────────────────────┘
                           ↓
┌──────────────────────────────────────────────────────────────────┐
│ STEP 6: Define Full Fit Function                                │
│                                                                   │
│   Template Fit Model:                                            │
│                                                                   │
│   C_HM(Δφ) = F × C_LM(Δφ)                                       │
│            + G × [1 + 2v₂²cos(2Δφ) + 2v₃²cos(3Δφ) + 2v₄²cos(4Δφ)]│
│                                                                   │
│   where:                                                         │
│   - C_HM: High multiplicity correlation (data)                  │
│   - C_LM: Low multiplicity correlation (template)               │
│   - F: Template normalization                                   │
│   - G: Flow amplitude                                            │
│   - v_n: Flow harmonics (what we want to extract!)             │
│                                                                   │
│   5 FREE PARAMETERS: [F, G, v₂², v₃², v₄²]                     │
└──────────────────────────────────────────────────────────────────┘
                           ↓
┌──────────────────────────────────────────────────────────────────┐
│ STEP 7: Perform RooFit Chi² Fit                                 │
│                                                                   │
│   Create RooFit objects:                                         │
│   ┌──────────────────────────────────────────────────────┐     │
│   │  RooRealVar F("F", "F", F_init, F_min, F_max);      │     │
│   │  RooRealVar G("G", "G", G_init, G_min, G_max);      │     │
│   │  RooRealVar v21("v21", "v21", 0.01, 0, 0.5);       │     │
│   │  RooRealVar v31("v31", "v31", 0.01, 0, 0.5);       │     │
│   │  RooRealVar v41("v41", "v41", 0.01, 0, 0.5);       │     │
│   │                                                       │     │
│   │  RooDataHist dsig("dsig", "dsig", x, hm);           │     │
│   │                                                       │     │
│   │  // Define PDF from template + flow function         │     │
│   │  RooAbsPdf* fitFunc = ...                            │     │
│   │                                                       │     │
│   │  // Perform chi² fit                                 │     │
│   │  fitFunc->chi2FitTo(dsig,                           │     │
│   │                     RooFit::SumW2Error(kTRUE),      │     │
│   │                     RooFit::PrintLevel(-1));        │     │
│   └──────────────────────────────────────────────────────┘     │
│                                                                   │
│   MINIMIZES:  χ² = Σᵢ [(yᵢ - fᵢ)² / σᵢ²]                      │
│                                                                   │
│   where σᵢ now includes systematics if LR_SYS_FRAC > 0!         │
└──────────────────────────────────────────────────────────────────┘
                           ↓
┌──────────────────────────────────────────────────────────────────┐
│ STEP 8: Calculate Chi²/ndf                                      │
│                                                                   │
│   RooAbsReal* chi2_obj = fitFunc->createChi2(dsig,             │
│                                               Extended(true),    │
│                                               SumW2);           │
│                                                                   │
│   chi2 = chi2_obj->getVal();                                    │
│   ndf = nBins - 5;                                              │
│   chi2ndf = chi2 / ndf;                                         │
│                                                                   │
│   Print:                                                         │
│   std::cout << "chi2: " << chi2 << std::endl;                  │
│   std::cout << "chi2/ndf: " << chi2ndf << std::endl;           │
└──────────────────────────────────────────────────────────────────┘
                           ↓
┌──────────────────────────────────────────────────────────────────┐
│ STEP 9: Quality Check                                           │
│                                                                   │
│   if (chi2ndf > 30000) {    // Template fit threshold           │
│       std::cout << "WARNING: Fit failed, chi2/ndf too large"   │
│                 << std::endl;                                   │
│                                                                   │
│       // Mark as failed:                                         │
│       F.setVal(-1);    F.setError(10);                          │
│       G.setVal(-1);    G.setError(10);                          │
│       v21.setVal(-1);  v21.setError(10);                        │
│       v31.setVal(-1);  v31.setError(10);                        │
│       v41.setVal(-1);  v41.setError(10);                        │
│       return FAILED;                                            │
│   }                                                              │
│                                                                   │
│   Good fits: 0.5 < χ²/ndf < 2                                   │
│   Acceptable: χ²/ndf < 10                                       │
│   Marginal: 10 < χ²/ndf < 30000                                │
│   Rejected: χ²/ndf ≥ 30000                                      │
└──────────────────────────────────────────────────────────────────┘
                           ↓
┌──────────────────────────────────────────────────────────────────┐
│ STEP 10: Extract Final Flow Coefficients                        │
│                                                                   │
│   Get fitted values:                                             │
│   v₂² = v21.getVal()    ±    σ_v₂² = v21.getError()            │
│   v₃² = v31.getVal()    ±    σ_v₃² = v31.getError()            │
│   v₄² = v41.getVal()    ±    σ_v₄² = v41.getError()            │
│                                                                   │
│   Convert to v_n:                                                │
│   v₂ = √(v₂²)           σ_v₂ = ½ × σ_v₂² / v₂                  │
│   v₃ = √(v₃²)           σ_v₃ = ½ × σ_v₃² / v₃                  │
│   v₄ = √(v₄²)           σ_v₄ = ½ × σ_v₄² / v₄                  │
│                                                                   │
│   NOTE: The errors σ_v_n include the effect of systematics      │
│         if LR_SYS_FRAC was set > 0 in Step 3!                  │
└──────────────────────────────────────────────────────────────────┘
                           ↓
                    ┌─────────────────┐
                    │ FINAL V_n       │
                    │ with errors     │
                    │ (stat + syst)   │
                    └─────────────────┘
```

---

## Systematic Uncertainty Addition

### Detailed Visualization

```
BEFORE Systematic Addition:          AFTER Systematic Addition:
(default, LR_SYS_FRAC=0.0)          (example: LR_SYS_FRAC=0.02, i.e., 2%)

Bin i:                               Bin i:
┌─────────────┐                      ┌─────────────┐
│ Content: N  │                      │ Content: N  │
│             │                      │             │
│  Error:     │                      │  Error:     │
│  σ = √N     │                      │  σ = √(N + (0.02×N)²) │
│             │                      │    = √(N + 0.0004N²)  │
│             │                      │    = √N × √(1 + 0.0004N) │
└─────────────┘                      └─────────────┘

Example values:                      Example with 2% systematic:
┌───────────┬─────────┬─────────┐   ┌───────────┬─────────┬─────────┬─────────┐
│ Bin │  N  │  σ_stat │         │   │ Bin │  N  │  σ_stat │  σ_sys  │ σ_total │
├───────────┼─────────┼─────────┤   ├───────────┼─────────┼─────────┼─────────┤
│  1  │  10 │   3.16  │         │   │  1  │  10 │   3.16  │  0.20   │  3.17   │
│  2  │ 100 │  10.0   │         │   │  2  │ 100 │  10.0   │  2.0    │  10.2   │
│  3  │1000 │  31.6   │         │   │  3  │1000 │  31.6   │  20.0   │  37.4   │
│  4  │5000 │  70.7   │         │   │  4  │5000 │  70.7   │ 100.0   │ 122.5   │
└───────────┴─────────┴─────────┘   └───────────┴─────────┴─────────┴─────────┘

Mathematical Formula:
σ_total = √(σ_stat² + σ_sys²)
        = √(N + (LR_SYS_FRAC × N)²)
        = √N × √(1 + LR_SYS_FRAC²× N)

For large N: σ_sys dominates
For small N: σ_stat dominates
```

### Impact on Chi² Calculation

```
WITHOUT systematics:                 WITH 2% systematics (LR_SYS_FRAC=0.02):

For bin with N=1000:                 For bin with N=1000:
                                    
data = 1000                          data = 1000
fit = 980                            fit = 980
σ = √1000 = 31.6                     σ = 37.4  (includes systematic)

χ²ᵢ = (1000-980)²/31.6²              χ²ᵢ = (1000-980)²/37.4²
    = 400/1000                           = 400/1399
    = 0.40                               = 0.29

Result: SMALLER contribution         Result: Systematic REDUCES the weight
        to total χ²                           of high-statistics bins

TOTAL χ² is SMALLER when              This can IMPROVE χ²/ndf for fits where
systematics are included              high-statistics bins dominate
```

---

## Numerical Examples

### Example 1: Simple Fourier Fit

**Input Data**: 50 bins, Δφ from -π/2 to 3π/2

| Bin | φ (center) | Data | Error | Fit | (Data-Fit)²/Error² |
|-----|-----------|------|-------|-----|---------------------|
| 1 | -1.57 | 100.5 | 10.0 | 98.2 | 0.0529 |
| 2 | -1.44 | 102.3 | 10.1 | 101.5 | 0.0063 |
| 3 | -1.32 | 104.1 | 10.2 | 105.1 | 0.0096 |
| ... | ... | ... | ... | ... | ... |
| 50 | 4.56 | 99.8 | 10.0 | 97.9 | 0.0361 |

**Sum**: χ² = 0.0529 + 0.0063 + 0.0096 + ... + 0.0361 = **48.3**

**Degrees of freedom**: ndf = 50 - 5 = **45**

**Reduced chi²**: χ²/ndf = 48.3/45 = **1.07** ✅ GOOD FIT

**Interpretation**: Data consistent with model within statistical fluctuations

---

### Example 2: Template Fit Without Systematics

**Configuration**: LR_SYS_FRAC = 0.0 (default)

| Parameter | Fitted Value | Error | Interpretation |
|-----------|--------------|-------|----------------|
| F | 0.987 | ±0.012 | Template normalization ≈ 1 |
| G | 52.3 | ±2.1 | Flow amplitude |
| v₂² | 0.0234 | ±0.0015 | → v₂ = 0.153 ±0.005 |
| v₃² | 0.0089 | ±0.0008 | → v₃ = 0.094 ±0.004 |
| v₄² | 0.0021 | ±0.0005 | → v₄ = 0.046 ±0.005 |

**Chi²/ndf**: 87.2 / 45 = **1.94** ✅ ACCEPTABLE

---

### Example 3: Template Fit With 1.5% Systematics

**Configuration**: LR_SYS_FRAC = 0.015

**Bin Error Comparison**:

| Bin | Content | σ_stat (old) | σ_sys | σ_total (new) | Change |
|-----|---------|--------------|-------|---------------|--------|
| 10 | 500 | 22.4 | 7.5 | 23.6 | +5.4% |
| 25 | 2000 | 44.7 | 30.0 | 53.9 | +20.6% |
| 40 | 5000 | 70.7 | 75.0 | 103.1 | +45.8% |

**Fit Results**:

| Parameter | Without Syst | With 1.5% Syst | % Change |
|-----------|--------------|----------------|----------|
| F | 0.987 ±0.012 | 0.986 ±0.014 | +16.7% error |
| G | 52.3 ±2.1 | 52.4 ±2.6 | +23.8% error |
| v₂² | 0.0234 ±0.0015 | 0.0235 ±0.0019 | +26.7% error |
| v₃² | 0.0089 ±0.0008 | 0.0089 ±0.0010 | +25.0% error |
| v₄² | 0.0021 ±0.0005 | 0.0021 ±0.0006 | +20.0% error |

**Chi²/ndf**: 87.2/45 = 1.94 → **61.3/45 = 1.36** ✅ IMPROVED

**Interpretation**: 
- Systematic uncertainties increase parameter errors by ~20-27%
- Chi²/ndf improves because high-statistics bins have less weight
- Central values remain similar (systematics don't bias the fit)

---

### Example 4: Failed Fit (Rejected)

**Scenario**: Poor template matching or wrong model

| Metric | Value | Status |
|--------|-------|--------|
| χ² | 1,453,782 | EXTREMELY HIGH |
| ndf | 45 | - |
| χ²/ndf | **32,306** | ❌ >> 30000 threshold |

**Action Taken**:
```cpp
if (chi2ndf > 30000) {
    // Set all parameters to -1 (failure flag)
    F = -1;    G = -1;
    v21 = -1;  v31 = -1;  v41 = -1;
    
    // Set errors to 10 (failure flag)
    return FAILED;
}
```

**Possible causes**:
- Template from wrong centrality bin
- Insufficient statistics
- Non-Gaussian fluctuations
- Wrong functional form
- Data corruption

---

## Summary Comparison Table

| Aspect | Fourier Fit | Template Fit (Default) | Template Fit (with Syst) |
|--------|-------------|------------------------|--------------------------|
| **Method** | Direct Fourier | Template subtraction | Template subtraction |
| **Parameters** | 5 (a₀...a₄) | 5 (F,G,v²ᵢ) | 5 (F,G,v²ᵢ) |
| **Systematics** | ❌ No | ❌ No (sysFrac=0) | ✅ Yes (sysFrac>0) |
| **Error Model** | σ = √N | σ = √N | σ = √(N+(frac×N)²) |
| **χ²/ndf Threshold** | 2000 | 30000 | 30000 |
| **Typical χ²/ndf** | 1-3 | 1-5 | 0.5-3 (often better) |
| **Use Case** | Simple correlations | Complex backgrounds | Publication-ready |
| **Runtime Config** | ❌ No | ✅ Yes (env vars) | ✅ Yes (env vars) |

---

## How to Enable Systematics

```bash
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# METHOD 1: No systematics (default)
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
root -l -b -q Process_TemplateFit.cxx

# Output: [FitDiag] NOT shown (no regularization applied)

#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# METHOD 2: Add 1% systematic
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
export LR_SYS_FRAC=0.01
root -l -b -q Process_TemplateFit.cxx

# Output: [FitDiag] Applying runtime error regularization: 
#         scale=1 sysFrac=0.01 floor=1e-06

#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# METHOD 3: Add 2% systematic + 10% error scaling
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
export LR_SYS_FRAC=0.02
export LR_ERROR_SCALE=1.1
root -l -b -q Process_TemplateFit.cxx

# Output: [FitDiag] Applying runtime error regularization: 
#         scale=1.1 sysFrac=0.02 floor=1e-06

#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# METHOD 4: Conservative - Large systematics
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
export LR_SYS_FRAC=0.05
export LR_ERROR_SCALE=1.2
export LR_ERROR_FLOOR=0.001
root -l -b -q Process_TemplateFit.cxx

# Output: [FitDiag] Applying runtime error regularization: 
#         scale=1.2 sysFrac=0.05 floor=0.001

#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Reset to defaults
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
unset LR_ERROR_SCALE
unset LR_SYS_FRAC
unset LR_ERROR_FLOOR
```

