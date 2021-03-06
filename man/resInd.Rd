% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resInd.R
\name{resInd}
\alias{resInd}
\title{Extracts change detection metrics from satellite time series}
\usage{
resInd(x, dates, type = "irregular", sc = 1, order = 3,
  formula = response ~ (trend + harmon), h = 0.15, plevel = 0.05, dr,
  drd, s = 3, NV = NA, plot = FALSE)
}
\arguments{
\item{x}{Numeric vector}

\item{dates}{Vector of dates in format "yyyy-mm-dd". Argument passed to bfastts() (see \code{\link[bfast]{bfastts}})}

\item{type}{Character. Time series type, "irregular", "16-day" or "10-day" (see \code{\link[bfast]{bfastts}}). Default="irregular"}

\item{sc}{Numeric. Scalefactor for input data. Default=1 (no rescaling). Scalefactor for input data. Default=1 (no rescaling).}

\item{order}{Numeric. Order of the harmonic component (see \code{\link[bfast]{bfastpp}}). Default=3.}

\item{formula}{Formula for the regression model and for fitting the breakpoints. Default=response ~ trend + harmon, a linear trend and a harmonic season component. Argument passed to breakpoints() and efp() (see \code{\link[strucchange]{breakpoints}} and  \code{\link[strucchange]{efp}}).}

\item{h}{Minimal fraction of oberservations required between two breakpoints relative to the sample size. Default=0.15. Argument passed to efp() and to breakpoints() (see \code{\link[strucchange]{breakpoints}} and  \code{\link[strucchange]{efp}}). Default=0.15}

\item{plevel}{Numeric value. Critical significance level for MOSUM test of structural stability. If test result is above the critical value, no breakpoints will be fitted. Default=0.05}

\item{dr}{Date Vector c(min, max) setting the time period during which a breakpoint is searched. Format: c(decimal year, decimal year); assuming 365 days/year}

\item{drd}{Date. Reference day (e.g.start of drought) based on which the 'timelag' is calculated. Format: decimal years (assuming 365 days/year)}

\item{s}{Numeric. Number of years used to calculate mean NDVI after the beginning of the time series used for calculation mean 'initial NDVI', and number of years before and after the breakpoint used to calculate 'PreNDVI', 'MagObsA' and 'MagObsR'. Assumes 365 days/year. Default=3}

\item{NV}{Numeric. Sets the value assigned to pixels without a breakpoint during the specified time window. Default: NA}

\item{plot}{Logical. If TRUE a plot is produced. Default: FALSE}
}
\value{
numeric vector ('resind') with 14 entries:
'BPNumb': Total number of breakpoints in time series
'Initial NDVI': Mean of data during the "s" first years of the time series.
'Intercept': Linear model intercept
'DBP': Drought Break Point yes/no (1/0). Yes, if a breakpoint occurs in time interval set with parameter "dr".
'BpTime': Timing of breakpoint. Format: Decimal year. If more than one breakpoints occurs during time interval, the first one is selected.
'Timelag': Number of days between drought reference day set with parameter "drd" and breakpoint
'RecTrend': Slope of linear trend in segment succeeding "drought breakpoint".
'PreTrend': Slope of linear trend in segment preceeding "drought breakpoint".
'PreNDVI': Mean of data of "s" years before drought breakpoint, based on observed data values.
'MagObsA': Absolute difference of mean observed data "s" years before and after the "drought breakpoint".
'MagObsR': Relative difference of mean observed data "s" years before and after the "drought breakpoint".
'MagTrendA': Absolute difference between last value of trend prediction before and first value of trend prediction after drought breakpoint. Based on corrected trend for irregular data. Still, with irregular data, the height of the trend line does not seem robust.
'MagTrendR': Relative difference between last value of trend prediction before and first value of trend prediction after drought breakpoint. Based on corrected trend for irregular data. Still, with irregular data, the height of the trend line does not seem robust.
'AmpDiffR': Relative difference in mean amplitudes (based on sine and cosine terms of harmonic model) in segment before and after drought breakpoint.
}
\description{
Computes several change detection metrics based on a BFAST type change detection framework for breakpoints ocurring during a given time period specified by the user (i.e. years of a known drought event). In addition, the overall number of breakpoints in the time series, the overall mean of the data and the inital value of the data at the beginning of the time series are extracted. Applicable to individual pixels.
}
\details{
This function was designed to explore several change detection metrics that can be extracted from a BFAST type change detection approach in relation to drought. The extracted metrics are at the experimental stage and should be used with caution. Not all metrics were found equally reliable: for an irregular time series the interecept of the linear trend line within segments (i.e. the height of the trend line after a breakpoint when plotted), was not stable. Therefore, the metrics relying directly on this information ('MagTrendA', 'MagTrendR') are equally not robust. The slope of the trend line within segments ('RecTrend','PreTrend') was found to be a robust model parameter, however, as well as the total number of breakpoints ('BPNumb').
Those parts of the code dealing with the fitting of breakpoints to an irregular time series, as well as the fitting of BFAST type models to a segmented time series was based on the function ("coefSegments" by Ben DeVries: \url{https://github.com/bendv/integrated-lts-cbm/blob/master/R/coefSegments.R}.
}
\examples{
#Load example raster data set (476 observations, 5x5 cells, including NA values) and date vector
data(stN) #raster brick
data(d) # date vector

#Plot first 9 layers of raster brick with NDVI scaling factor
sc <- 0.0001
library(raster)
plot(stN*sc, 1:9)

##With package "bfastSpatial" you can extract information on acquisition date
# from typical Landsat file names \\url{https://github.com/loicdtx/bfastSpatial}
#gs <- getSceneinfo(names(stN))
#d <- gs$date

#Select target pixel from raster stack
targcell <- 1

#Extract vector of NDVI values for targcell and plot time series.
x <- as.vector(stN[targcell])
plot(d,x*sc, ylab='NDVI', xlab='year', main='targcell', ylim=c(0,1))

#Run function "resInd"
y <- resInd(x, d, dr=c(2004.753,2008.751), drd=2004.753, plot=TRUE)

#Should return a plot and a vector y containing 14 values


}
\author{
Jennifer von Keyserlingk
}
