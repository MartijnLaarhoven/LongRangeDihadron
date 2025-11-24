opt="-b --configuration json://configuration_longRangeDihadronCor.json"
    aodfile="AO2D_LHC25ae_apass1_564400.root"   
    o2-analysis-propagationservice ${opt} |
    o2-analysis-trackselection ${opt} |
    o2-analysis-pid-tof-full ${opt} |
    o2-analysis-pid-tof-base ${opt} |
    o2-analysis-pid-tof-beta ${opt} |
    o2-analysis-pid-tpc-service ${opt} |
    o2-analysis-pid-its ${opt} |
    o2-analysis-multcenttable ${opt} |
    o2-analysis-event-selection-service ${opt} |
    o2-analysis-ft0-corrected-table ${opt} |
    o2-analysis-cf-long-range-dihadron-cor ${opt}  --aod-file ${aodfile}
