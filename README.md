# resInd
Code for extraction of change detection metrics based on BFAST framework.
Please note that the extracted change detection metrics are at an experimental stage and were designed for exploration porposes.

For the manuscript by von Keyserlingk and de Hoop et al. (currently under revision), we have applied the function resIndSpatial.R to our study area. For analysis, we have used the following parameters:

BPNumb: to extract the total number of breakpoints fitted to the NDVI time series. We used this as a measure for vegetation's long-term resistance to climatic variations

DBP: to check if a breakpoint occurred during the hydrological years 2005 & 2008. 

MagObsR: if a breakpoint was found during this time, we used MagObsR to extract the relative difference in NDVI before and after this breakpoint.

RecTrend: If MagOBSR was below or equal to -10%, we used RecTrend to estimate the linear "NDVI recovery trend" after the breakpoint. We used this as a measure for the recovery rate after drought.

PreNDVI: This was used as a measure for the NDVI before the breakpoint (used in the regression analysis)



