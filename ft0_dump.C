// Helper macro: dump FT0Amp / FT0AmpCorrect per-channel entries and mean amplitudes.
// Usage:
//   root -l -q 'ft0_dump.C("../AnalysisResultsROOTFiles/longRangeDihadronCorr/AnalysisResults_LHC25af_pass2_565246.root", "long-range-dihadron-cor_Cent_0_20")'
// Optional:
//   third arg = FT0Amp histogram name (default "FT0Amp")
//   fourth arg = FT0AmpCorrect histogram name (default "FT0AmpCorrect")
// Notes:
//   Channels 0-95 = FT0A, 96-207 = FT0C (as produced in longRangeDihadronCor).
//   The macro prints: channel, side, entries, mean amplitude, mean corrected amplitude (if available).

#include <TFile.h>
#include <TH2.h>
#include <TProfile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>

void ft0_dump(const char* filename = "AnalysisResults.root",
              const char* taskDir = "long-range-dihadron-cor",
              const char* hName = "FT0Amp",
              const char* hNameCorr = "FT0AmpCorrect",
              const char* pdfOut = "")
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

  // Projections and profiles
  auto hEntries = h->ProjectionX("ft0_ch_entries");      // entries per channel
  auto hMean    = h->ProfileX("ft0_ch_mean");           // mean amplitude per channel
  TProfile* hMeanCorr = hCorr ? hCorr->ProfileX("ft0_ch_mean_corr") : nullptr;
  std::vector<std::string> blockSummaries;

  // The histogram is booked with 220 bins (0..219) while FT0 channels are 0..207.
  // Restrict the loop to physical channels unless higher bins carry content.
  const int nBins = hEntries->GetNbinsX();
  const int maxChIdx = 207; // highest real FT0 channel id
  const int loopBins = std::min(nBins, maxChIdx + 1);
  std::cout << std::left << std::setw(6) << "ch" << std::setw(6) << "side"
            << std::setw(12) << "entries" << std::setw(16) << "meanAmp"
            << std::setw(16) << "meanAmpCorr" << "\n";
  for (int i = 1; i <= loopBins; i += 4) {
    double blockEntries = 0.0;
    double blockMeanAmp = 0.0;
    double blockMeanAmpCorr = 0.0;
    int blockCount = 0;
    for (int j = 0; j < 4 && (i + j) <= loopBins; ++j) {
      int bin = i + j;
      int ch = bin - 1; // channel ID matches bin center index offset by 1
      const char* side = (ch < 96) ? "FT0A" : "FT0C";
      double entries = hEntries->GetBinContent(bin);
      double mean = hMean ? hMean->GetBinContent(bin) : 0.0;
      double meanCorr = hMeanCorr ? hMeanCorr->GetBinContent(bin) : 0.0;
      blockEntries += entries;
      blockMeanAmp += mean * entries; // weight by entries to combine sensibly
      if (hMeanCorr)
        blockMeanAmpCorr += meanCorr * entries;
      blockCount++;
      std::cout << std::left << std::setw(6) << ch
                << std::setw(6) << side
                << std::setw(12) << entries
                << std::setw(16) << mean
                << std::setw(16) << (hMeanCorr ? meanCorr : NAN)
                << "\n";
    }
    // Print combined stats for the 4-channel block (entries-weighted mean)
    double avgMean = (blockEntries > 0) ? (blockMeanAmp / blockEntries) : 0.0;
    double avgMeanCorr = (hMeanCorr && blockEntries > 0) ? (blockMeanAmpCorr / blockEntries) : NAN;
    int firstCh = i - 1;
    int lastCh = std::min(i + 3, loopBins) - 1;
    std::cout << "-- block " << firstCh << "-" << lastCh
              << " entries=" << blockEntries
              << " meanAmp=" << avgMean
              << " meanAmpCorr=" << (hMeanCorr ? avgMeanCorr : NAN)
              << "\n\n";
    std::cout << std::flush;

    std::ostringstream os;
    os << "block " << firstCh << "-" << lastCh
       << " entries=" << blockEntries
       << " meanAmp=" << avgMean
       << " meanAmpCorr=" << (hMeanCorr ? avgMeanCorr : NAN);
    blockSummaries.push_back(os.str());
  }

  // Optional PDF output of block summaries (paginated)
  if (pdfOut && std::string(pdfOut).size() > 0) {
    const int linesPerPage = 32;
    int total = static_cast<int>(blockSummaries.size());
    TCanvas c("c_ft0_dump", "ft0_dump", 900, 1200);
    for (int start = 0; start < total; start += linesPerPage) {
      c.Clear();
      TPaveText pt(0.05, 0.05, 0.95, 0.95, "NDC");
      pt.SetTextFont(42);
      pt.SetTextSize(0.03);
      pt.SetBorderSize(0);
      pt.SetFillStyle(0);
      for (int k = 0; k < linesPerPage && (start + k) < total; ++k) {
        pt.AddText(blockSummaries[start + k].c_str());
      }
      pt.Draw();
      if (start == 0 && start + linesPerPage < total) c.Print((std::string(pdfOut) + "(").c_str());
      else if (start + linesPerPage >= total && start == 0) c.Print(pdfOut);
      else if (start + linesPerPage >= total) c.Print((std::string(pdfOut) + ")").c_str());
      else c.Print(pdfOut);
    }
  }

  f->Close();
}
