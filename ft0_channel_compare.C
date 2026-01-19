// Enhanced FT0 channel analysis macro: produces a table comparing each channel to its mirror, with improved formatting and summary statistics.
// Usage:
//   root -l -q 'ft0_channel_compare.C("../AnalysisResultsROOTFiles/longRangeDihadronCorr/AnalysisResults_LHC25af_pass2_565246.root", "long-range-dihadron-cor_Cent_0_20", "ft0_channel_compare.csv")'
//
// Output: CSV file with columns: ch,side,entries,meanAmp,meanAmpCorr,mirror_ch,mirror_side,mirror_entries,mirror_meanAmp,mirror_meanAmpCorr,diff_meanAmp,diff_meanAmpCorr

#include <TFile.h>
#include <TH2.h>
#include <TProfile.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

void ft0_channel_compare(const char* filename = "AnalysisResults.root",
                        const char* taskDir = "long-range-dihadron-cor",
                        const char* csvOut = "ft0_channel_compare.csv",
                        const char* hName = "FT0Amp",
                        const char* hNameCorr = "FT0AmpCorrect")
{
  auto f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[Error] Cannot open file: " << (filename ? filename : "(null)") << std::endl;
    return;
  }

  std::string path = std::string(taskDir) + "/" + hName;
  auto h = dynamic_cast<TH2*>(f->Get(path.c_str()));
  if (!h) {
    std::cerr << "[Error] Histogram not found: " << path << std::endl;
    f->Close();
    return;
  }

  std::string pathCorr = std::string(taskDir) + "/" + hNameCorr;
  auto hCorr = dynamic_cast<TH2*>(f->Get(pathCorr.c_str()));

  auto hEntries = h->ProjectionX("ft0_ch_entries");
  auto hMean    = h->ProfileX("ft0_ch_mean");
  TProfile* hMeanCorr = hCorr ? hCorr->ProfileX("ft0_ch_mean_corr") : nullptr;

  const int nCh = 208; // 0..207
  std::ofstream csv(csvOut);
  csv << "ch,side,entries,meanAmp,meanAmpCorr,mirror_ch,mirror_side,mirror_entries,mirror_meanAmp,mirror_meanAmpCorr,diff_meanAmp,diff_meanAmpCorr\n";

  for (int ch = 0; ch < nCh; ++ch) {
    int mirror_ch = (ch < 96) ? (207 - ch) : (207 - ch); // FT0A: 0-95, FT0C: 96-207
    const char* side = (ch < 96) ? "FT0A" : "FT0C";
    const char* mirror_side = (mirror_ch < 96) ? "FT0A" : "FT0C";
    int bin = ch + 1;
    int mirror_bin = mirror_ch + 1;
    double entries = hEntries->GetBinContent(bin);
    double mean = hMean->GetBinContent(bin);
    double meanCorr = hMeanCorr ? hMeanCorr->GetBinContent(bin) : 0.0;
    double mirror_entries = hEntries->GetBinContent(mirror_bin);
    double mirror_mean = hMean->GetBinContent(mirror_bin);
    double mirror_meanCorr = hMeanCorr ? hMeanCorr->GetBinContent(mirror_bin) : 0.0;
    double diff_mean = mean - mirror_mean;
    double diff_meanCorr = meanCorr - mirror_meanCorr;
    csv << ch << "," << side << "," << entries << "," << mean << "," << meanCorr << ","
        << mirror_ch << "," << mirror_side << "," << mirror_entries << "," << mirror_mean << "," << mirror_meanCorr << ","
        << diff_mean << "," << diff_meanCorr << "\n";
  }
  csv.close();
  std::cout << "[Info] Channel comparison table written to: " << csvOut << std::endl;
  f->Close();
}
