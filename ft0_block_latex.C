// Macro to produce a LaTeX table of FT0 block averages, showing each block and its mirrored block side-by-side.
// Usage:
//   root -l -q 'ft0_block_latex.C("../AnalysisResultsROOTFiles/longRangeDihadronCorr/AnalysisResults_LHC25ae_pass2_573730.root", "long-range-dihadron-cor", "ft0_block_latex.tex")'
// Output: LaTeX file with block averages and mirrored block comparison.

#include <TFile.h>
#include <TH2.h>
#include <TProfile.h>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

// Usage:
//   root -l -q 'ft0_block_latex.C("../AnalysisResultsROOTFiles/longRangeDihadronCorr/AnalysisResults_LHC25ae_pass2_573730.root", "long-range-dihadron-cor_Cent_0_20", "ft0_block_latex.tex", "FT0Amp", "FT0AmpCorrect")'
//   or specify histogram names:
//   root -l -q 'ft0_block_latex.C("../AnalysisResultsROOTFiles/longRangeDihadronCorr/AnalysisResults_LHC25ae_pass2_573730.root", "long-range-dihadron-cor", "ft0_block_latex.tex", "FT0Amp", "FT0AmpCorrect")'
void ft0_block_latex(const char* filename = "AnalysisResults.root",
                    const char* taskDir = "long-range-dihadron-cor",
                    const char* texOut = "ft0_block_latex.tex",
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
  if (!hCorr) {
    std::cerr << "[Info] Corrected histogram not found (will skip): " << pathCorr << std::endl;
  }

  auto hEntries = h->ProjectionX("ft0_ch_entries");
  auto hMean    = h->ProfileX("ft0_ch_mean");
  TProfile* hMeanCorr = hCorr ? hCorr->ProfileX("ft0_ch_mean_corr") : nullptr;

  const int nCh = 208; // 0..207
  const int blockSize = 4;
  std::ofstream tex(texOut);
  tex << "\\begin{table}[htbp]\n";
  tex << "\\centering\n";
  tex << "\\caption{FT0 block averages and mirrored block comparison}" << std::endl;
  tex << "\\label{tab:ft0blocks}" << std::endl;
  tex << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}\\hline\n";
  tex << "Block & Side & Entries & MeanAmp & MeanAmpCorr & Mirrored Block & Side & Entries & MeanAmp & MeanAmpCorr \\" << std::endl;
  tex << "\\hline" << std::endl;

  for (int block = 0; block < nCh; block += blockSize) {
    // Regular block
    double blockEntries = 0, blockMeanAmp = 0, blockMeanAmpCorr = 0;
    int blockFirstCh = block, blockLastCh = std::min(block + blockSize - 1, nCh - 1);
    for (int ch = blockFirstCh; ch <= blockLastCh; ++ch) {
      int bin = ch + 1;
      double entries = hEntries->GetBinContent(bin);
      double mean = hMean->GetBinContent(bin);
      double meanCorr = hMeanCorr ? hMeanCorr->GetBinContent(bin) : 0.0;
      blockEntries += entries;
      blockMeanAmp += mean * entries;
      if (hMeanCorr) blockMeanAmpCorr += meanCorr * entries;
    }
    double avgMean = (blockEntries > 0) ? (blockMeanAmp / blockEntries) : 0.0;
    double avgMeanCorr = (hMeanCorr && blockEntries > 0) ? (blockMeanAmpCorr / blockEntries) : 0.0;
    std::string side = (blockFirstCh < 96) ? "FT0A" : "FT0C";

    // Mirrored block
    int mirrorBlockFirstCh = nCh - blockLastCh - 1;
    int mirrorBlockLastCh = nCh - blockFirstCh - 1;
    double mirrorBlockEntries = 0, mirrorBlockMeanAmp = 0, mirrorBlockMeanAmpCorr = 0;
    for (int ch = mirrorBlockFirstCh; ch <= mirrorBlockLastCh; ++ch) {
      int bin = ch + 1;
      double entries = hEntries->GetBinContent(bin);
      double mean = hMean->GetBinContent(bin);
      double meanCorr = hMeanCorr ? hMeanCorr->GetBinContent(bin) : 0.0;
      mirrorBlockEntries += entries;
      mirrorBlockMeanAmp += mean * entries;
      if (hMeanCorr) mirrorBlockMeanAmpCorr += meanCorr * entries;
    }
    double mirrorAvgMean = (mirrorBlockEntries > 0) ? (mirrorBlockMeanAmp / mirrorBlockEntries) : 0.0;
    double mirrorAvgMeanCorr = (hMeanCorr && mirrorBlockEntries > 0) ? (mirrorBlockMeanAmpCorr / mirrorBlockEntries) : 0.0;
    std::string mirrorSide = (mirrorBlockFirstCh < 96) ? "FT0A" : "FT0C";

    // Escape underscores in side and block labels if present
    auto escape = [](const std::string& s) {
      std::string out;
      for (char c : s) {
        if (c == '_') out += "\\_";
        else out += c;
      }
      return out;
    };
    tex << escape(std::to_string(blockFirstCh) + "-" + std::to_string(blockLastCh)) << " & " << escape(side) << " & " << blockEntries << " & " << avgMean << " & " << avgMeanCorr << " & "
      << escape(std::to_string(mirrorBlockFirstCh) + "-" + std::to_string(mirrorBlockLastCh)) << " & " << escape(mirrorSide) << " & " << mirrorBlockEntries << " & " << mirrorAvgMean << " & " << mirrorAvgMeanCorr << " \\" << std::endl;
    tex << "\\hline" << std::endl;
  }
  tex << "\\end{tabular}\n";
  tex << "\\end{table}\n";
  tex.close();
  std::cout << "[Info] LaTeX block table written to: " << texOut << std::endl;
  f->Close();
}
