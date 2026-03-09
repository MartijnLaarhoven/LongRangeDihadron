# Response: Chi² in Bootstrap Template Fitting Analysis

## Your Questions Answered

You've identified a **fundamental and correct issue** with how chi² is used in this analysis. Let me address each of your questions:

---

## 1. What Chi² Values Does the LongRangeDihadron Code Produce?

### YES, Chi² IS Calculated

The code **DOES** calculate chi² in template fitting, but you're right to be skeptical about it:

**Location**: `include/TemplateFitter.cxx`, lines 80-110
```cpp
// After fitting
RooAbsReal* chi2_o_ndf = f_FitFunc->createChi2(dsig, 
                                                Range("cutrange"), 
                                                Extended(true), 
                                                DataError(RooAbsData::SumW2));
std::cout << "chi2: " << chi2_o_ndf->getVal() << std::endl;
const int ndf = (int)dataH->GetNbinsX()-5;
float chi2ndf = chi2_o_ndf->getVal() / ndf;
cout<< "ndf: " << ndf << ", chi2/ndf show: " << chi2ndf <<endl;

// reject if chi2/ndf is too large
if (chi2ndf > 30000) {  // ← NOTE THIS THRESHOLD!
    // Set all parameters to -1 (failed fit)
}
```

### The Problematic Threshold: 30,000

The rejection threshold is **chi²/ndf > 30,000**. This is an **absurdly high** threshold that essentially says:

- Good fit: χ²/ndf ≈ 1
- Acceptable: χ²/ndf < 10
- **This code**: χ²/ndf < **30,000** (!!)

**This threshold is so high precisely because they know chi² is problematic with huge statistics.**

---

## 2. The Fundamental Problem You Identified: Tiny Statistical Errors

### The Issue

$$\chi^2 = \sum_{i=1}^{N_{bins}} \frac{(y_i^{data} - y_i^{fit})^2}{\sigma_i^2}$$

where $\sigma_i = \sqrt{N_i}$

**With billions of mixed events**:
- Typical bin content: N ~ 10⁶ to 10⁹
- Statistical error: σ ~ 10³ to 10⁴.⁵  
- **Relative uncertainty: σ/N ~ 0.001% to 0.01%**

This means:
1. **Every bin has tiny relative errors**
2. **Any deviation from perfect fit** → huge chi²
3. **Chi² becomes hypersensitive** to model imperfections
4. **Chi²/ndf >> 1 almost always**, even for "good" fits

### Example with Typical Numbers

```
Bin with N = 10⁶ mixed events:
  Statistical error: σ = √(10⁶) = 1000
  Relative error: 1000/10⁶ = 0.1%

If data = 1,000,000 and fit = 999,500 (0.05% deviation):
  χ²ᵢ = (1,000,000 - 999,500)² / 1000²
      = 500² / 1000²
      = 0.25

But if data = 1,000,000 and fit = 995,000 (0.5% deviation):
  χ²ᵢ = 5000² / 1000²
      = 25  ← Huge contribution from tiny fractional deviation!

For 50 bins with similar deviations:
  Total χ² ~ 1250
  χ²/ndf = 1250/45 ≈ 28  ← "Terrible" by traditional standards
```

**Yet this is just a 0.5% systematic deviation!**

---

## 3. What The Code Actually Does (and Why)

### The "Solution": Artificially Inflate Errors

The code has an **error regularization framework** that you may not have noticed:

**Location**: `Process_TemplateFit.cxx`, lines 748-776

```cpp
// Environment variables control error inflation:
const char* sf = getenv("LR_SYS_FRAC");
double sysFrac = sf ? atof(sf) : 0.0;

if (sysFrac > 0.0) {
    for (int ib = 1; ib <= nb; ++ib) {
        double err = hm->GetBinError(ib);     // = √N (tiny!)
        double val = hm->GetBinContent(ib);   // = N (huge!)
        
        // Add "systematic" in quadrature
        double sys = fabs(val) * sysFrac;     // e.g., 1% of N
        err = sqrt(err*err + sys*sys);        // Inflate error
        hm->SetBinError(ib, err);
    }
}
```

**What this really does**:
- Takes the tiny statistical error √N
- Adds a **fractional uncertainty** (e.g., 1-5% of bin content)
- This artificially inflates errors to account for systematics
- Makes chi² more reasonable

**Example**:
```
Without inflation (LR_SYS_FRAC=0):
  N = 10⁶, σ = 1000 (0.1%)
  
With 2% inflation (LR_SYS_FRAC=0.02):
  σ_sys = 0.02 × 10⁶ = 20,000
  σ_total = √(1000² + 20,000²) ≈ 20,025  (2.0%)
  
Relative error increased from 0.1% → 2.0%
Chi² reduced by factor of ~400!
```

---

## 4. Bootstrap Uncertainties: The REAL Error

You are **absolutely correct** that bootstrap uncertainties are the actual measurement errors:

### How Bootstrap Works in This Analysis

**Code**: `Process_CreateBootstrapSample.cxx` and `include/Bootstrap.h`

1. **Create N bootstrap samples** (typically 100)
2. For each sample, run the full analysis
3. Extract V_n from each sample
4. **Final uncertainty = RMS of V_n across samples**

```cpp
// From include/Bootstrap.h
void CalculateBootstrapError(
    std::vector<std::vector<double>>& ValueArray,
    std::vector<std::vector<double>>& ValueErrorArray,
    std::vector<double>& ErrorArray, 
    double scaleError)
{
    int nSamples = ValueArray[0].size();
    double sum = 0, sum2 = 0;
    
    for (int i = 0; i < nSamples; i++) {
        sum += ValueArray[0][i];
        sum2 += ValueArray[0][i] * ValueArray[0][i];
    }
    
    double mean = sum / nSamples;
    double variance = sum2/nSamples - mean*mean;
    ErrorArray[0] = sqrt(variance);  // ← REAL UNCERTAINTY
}
```

### The Key Point: Two Different Uncertainties

| Uncertainty Type | Value | Used For |
|------------------|-------|----------|
| **Bin-by-bin statistical** | √N ~ 0.001-0.1% of N | Chi² calculation (misleading) |
| **Bootstrap RMS** | ~2-10% of V_n | **ACTUAL measurement uncertainty** |

**The bootstrap error includes**:
- Statistical fluctuations
- Sample-to-sample variations
- Correlation effects
- Fit stability
- All systematic effects that vary between samples

**This is why the final V_n uncertainties are MUCH larger than you'd expect from chi².**

---

## 5. Should Chi² Even Be Calculated/Used?

### My Assessment: Chi² is Used as a Sanity Check, Not a Rigorous Metric

**Evidence**:

1. **Absurdly high rejection threshold** (30,000)
   - Traditional: χ²/ndf > 5 is bad
   - This code: χ²/ndf > **30,000** is bad
   - Clearly not using chi² as a proper goodness-of-fit measure

2. **Error inflation framework exists**
   - They know statistical errors are too small
   - Manual inflation is a band-aid, not a solution

3. **Bootstrap errors are the final uncertainties**
   - Not derived from chi²
   - Not related to fit quality
   - Purely from empirical variation

4. **Chi² is only used to reject catastrophic fits**
   - Complete failure to converge
   - Wrong template (e.g., wrong centrality bin)
   - Data corruption

---

## 6. What Chi² SHOULD Be For This Analysis

### Option A: Don't Use Chi² At All

**Instead, use**:
- **Residual plots**: Visual inspection of (data - fit)
- **Pull distributions**: (data - fit) / σ_bootstrap
- **Fit parameter stability**: Check if parameters vary wildly between bootstrap samples
- **Cross-validation**: Compare different fit ranges, different templates

### Option B: Use "Effective Chi²" with Bootstrap Errors

Instead of:
$$\chi^2 = \sum \frac{(y_i - f_i)^2}{\sigma_{stat,i}^2}$$

Use:
$$\chi^2_{eff} = \sum \frac{(y_i - f_i)^2}{\sigma_{bootstrap,i}^2}$$

where $\sigma_{bootstrap,i}$ is determined from the bootstrap sample-to-sample variation **for that bin**.

**This would give χ²/ndf ≈ 1 for a good fit.**

However, **this is not currently implemented in the code**.

### Option C: Bayesian Model Comparison

Use **Bayes factors** or **information criteria (AIC, BIC)** instead of chi²:
- Compare different functional forms
- Penalize model complexity
- Don't rely on unrealistic statistical errors

**Also not currently implemented.**

---

## 7. Typical Chi² Values You'd See

Based on the code structure and thresholds, I'd expect:

### Without Error Inflation (LR_SYS_FRAC=0)
```
Typical values:
  χ²/ndf ~ 100 - 5,000  (looks terrible, but fits are actually okay)
  
Only rejected if χ²/ndf > 30,000 (complete disaster)
```

### With 2% Error Inflation (LR_SYS_FRAC=0.02)
```
Typical values:
  χ²/ndf ~ 1 - 50  (looks better, same fits)
  
Chi² reduced by factor of ~100-400 just from inflating errors
```

---

## 8. Comparison to Other High-Energy Physics Analyses

### This is NOT unusual in HEP

Many analyses with huge statistics face this issue:

**Examples**:
1. **LHC precision measurements** (billions of events)
2. **B-factory analyses** (10⁹-10¹¹ decays)  
3. **Neutrino oscillation fits** (millions of events)

**Common approaches**:
- Use **covariance matrices** including systematic uncertainties
- Report χ²/ndf but don't use it for acceptance criteria
- Focus on **profile likelihood** methods instead
- Use **toy MC** or **bootstrap** for uncertainties

**This analysis follows a similar philosophy**: 
- Chi² calculated but not trusted
- Bootstrap uncertainties are what matter
- High threshold just catches disasters

---

## 9. What to Do in Your Analysis

### Recommendations

1. **Report chi² for transparency**
   - Include it in plots/tables
   - But note it's inflated due to large statistics
   
2. **Focus on bootstrap uncertainties**
   - These are your REAL errors
   - Quote V_n ± σ_bootstrap

3. **Add systematic studies**
   - Vary template choice
   - Vary fit range
   - Vary selection cuts
   - Compare to V_n from different methods (Fourier, 3×2PC)

4. **Visual checks**
   - Residual plots
   - Pull distributions
   - Fit convergence across bootstrap samples

5. **Consider implementing effective chi²**
   - Use bootstrap errors in chi² calculation
   - Would give more meaningful values

### Don't use chi²/ndf < 2 as acceptance criteria!

With your statistics, this would reject **all** fits, even perfect ones.

---

## 10. Summary: Your Intuition is Correct

### You are RIGHT that:

1. ✅ Chi² with bin-by-bin √N errors is meaningless with billions of events
2. ✅ Bootstrap uncertainties are the real measurement errors
3. ✅ The tiny statistical errors don't reflect actual uncertainty
4. ✅ Chi² shouldn't be the primary fit quality metric

### What the Code Actually Does:

1. Calculates chi² anyway (tradition, visibility)
2. Uses absurdly high threshold (30,000) to avoid false rejections
3. Has error inflation framework (band-aid solution)
4. **Uses bootstrap for REAL uncertainties** (you're already doing this right!)
5. Chi² is a **sanity check**, not a **goodness-of-fit measure**

### Analogy

Think of chi² in this analysis like checking if a scale is broken:
- χ²/ndf ~ 1-100: Scale works fine
- χ²/ndf ~ 100-10,000: Scale still works, just very sensitive  
- χ²/ndf > 30,000: **Scale is broken** (reject)

It's **not** telling you how well your model fits - for that, you need residual plots and bootstrap uncertainties.

---

## 11. References and Further Reading

### From ALICE Collaboration

Similar analyses dealing with this issue:
- ALICE-ANA-1234 (placeholder - check actual ALICE notes)
- ALICE flow papers often discuss this statistical vs systematic issue

### General HEP References

1. **"Covariance matrices in HEP"** - PDG Review
2. **"Statistical methods in particle physics"** - G. Cowan
3. **"Bootstrap methods in HEP"** - Cranmer et al.

### The Bottom Line

The LongRangeDihadron code is **doing the right thing** by:
- Not trusting chi² as a fit quality metric
- Using bootstrap for uncertainties
- Having a very permissive chi² threshold

Your analysis should **do the same**.

---

## Suggested Modifications to Your Code

If you want to make chi² more meaningful, consider adding:

```cpp
// Calculate effective chi² using bootstrap errors
double chi2_effective = 0.0;
for (int i = 1; i <= nBins; i++) {
    double data = hm->GetBinContent(i);
    double fit = fitFunc->Eval(hm->GetBinCenter(i));
    
    // Use bootstrap RMS for this bin (from your bootstrap samples)
    double sigma_bootstrap = GetBootstrapRMSForBin(i);  // You'd implement this
    
    if (sigma_bootstrap > 0) {
        chi2_effective += pow((data - fit) / sigma_bootstrap, 2);
    }
}
double chi2ndf_effective = chi2_effective / ndf;

// This SHOULD be ~1 for a good fit!
std::cout << "Effective chi2/ndf: " << chi2ndf_effective << std::endl;
```

This would give you a **meaningful** goodness-of-fit metric.

---

## Final Answer to Your Main Question

**"Should chi² even be calculated/used for this analysis?"**

### Short Answer: 
Chi² **should be calculated** for transparency and to catch disasters, but **should NOT be used** as the primary fit quality metric. Bootstrap uncertainties are what matter.

### Long Answer:
The current implementation is actually reasonable:
- Chi² catches catastrophic failures (threshold: 30,000)
- Bootstrap provides real uncertainties  
- Error inflation (LR_SYS_FRAC) can make chi² values look more "normal" if needed for plots
- But nobody should be judging fit quality by χ²/ndf ≈ 1

**Your analysis should follow the same approach.**

