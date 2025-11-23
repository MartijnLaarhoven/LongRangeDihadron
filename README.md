# DihadronAnalysis

## Description
The basic post-processing framework for TPC-FT0A/C di-hadron correlations, written by Zhiyong (zhiyong.lu@cern.ch).

This framework works for the output from the long-range-dihadron-cor task (O2Physics/PWGCF/TwoParticleCorrelations/Tasks/longRangeDiHadronCor.cxx)

## Usage
Should be run in the following sequence to get the results:

- root -l -b -q Process_dPhidEta.cxx
- root -l -b -q Process_CreateBootstrapSample.cxx
- root -l -b -q Process_TemplateFit.cxx

If you need 3 times 2PC results:
- root -l -b -q Process_3times2PC.cxx
# LongRangeDihadron
