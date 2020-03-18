# resInd
Code for extraction of change detection metrics based on BFAST framework.
Please note that some of the change detection metrics returned by the functions are at an experimental stage and were designed for exploration purposes.

For the manuscript by von Keyserlingk and de Hoop et al. (currently under revision), we have applied the function "resIndSpatial.R" to our study area. For this analysis, we have only used the following parameters, which we consider robust:

BPNumb: to extract the total number of breakpoints fitted to the NDVI time series. We used this as a measure for vegetation's long-term resistance to climatic variations.

DBP: to check if a breakpoint occurred between the hydrological years 2005 & 2008. 

MagObsR: if DBP was true, we used MagObsR to extract the relative difference in NDVI before and after this breakpoint.

RecTrend: If MagOBSR was below or equal to -10%, we used RecTrend to estimate the linear "NDVI recovery trend" after the breakpoint. We used this as a measure for the recovery rate after drought.

PreNDVI: This was used as a measure for the "NDVI before breakpoint" (used in the regression analysis).



Installation:
The package can be installed from github with the package "devtools":

library(devtools)
install_github('jennifervk/resInd')



