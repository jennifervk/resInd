
#' @title Runs function "resInd" on gridded time series and returns spatial output
#' @description Computes several change detection metrics based on a BFAST (Break for Additive Season and Trend)
#' change detection framework. These change detection metrics are calculated for breakpoints occurring during a
#' given time period specified by the user (i.e. during a known drought). In addition, data about the full time
#' series is extracted: the overall number of breakpoints as well as the overall mean and initial value of the
#' dependent time series data. The included functions were designed for exploring the use of a breakpoint model
#' to study the response of vegetation to drought and climatic perturbartions. Please note that some of the
#' extracted change detection metrics are at an experimental stage. Applicable to to gridded time-series data.
#' @param x RasterBrick or rasterStack object
#' @param dates See \code{\link{resInd}}
#' @param type See \code{\link{resInd}}
#' @param sc See \code{\link{resInd}}
#' @param order See \code{\link{resInd}}
#' @param formula See \code{\link{resInd}}
#' @param h See \code{\link{resInd}}
#' @param plevel See \code{\link{resInd}}
#' @param dr See \code{\link{resInd}}
#' @param drd See \code{\link{resInd}}
#' @param s See \code{\link{resInd}}
#' @param NV See \code{\link{resInd}}
#' @param mc.cores Numeric. Number of cores used for the calcualtion. Default=1
#'
#' @return Large Raster stack with 14 layers:
#' 'BPNumb': Total number of breakpoints in time series
#' 'Initial NDVI': Mean of data during the "s" first years of the time series.
#' 'Intercept': Linear model intercept
#' 'DBP': Drought Break Point yes/no (1/0). Yes, if a breakpoint occurs in time
#'   interval set with parameter "dr".
#' 'BpTime': Timing of breakpoint. Format: Decimal year. If more than one breakpoints
#'   occurs during time interval, the first one is selected.
#' 'Timelag': Number of days between drought reference day set with parameter "drd"
#'   and breakpoint
#' 'RecTrend': Slope of linear trend in segment succeeding "drought breakpoint".
#' 'PreTrend': Slope of linear trend in segment preceeding "drought breakpoint".
#' 'PreNDVI': Mean of data of "s" years before drought breakpoint, based on observed
#'  data values.
#' 'MagObsA': Absolute difference of mean observed data "s" years before and after
#'   the "drought breakpoint".
#' 'MagObsR': Relative difference of mean observed data "s" years before and after
#'   the "drought breakpoint".
#' 'MagTrendA': Absolute difference between last value of trend prediction before
#'   and first value of trend prediction after drought breakpoint. Based on corrected
#'   trend for irregular data. Still, with irregular data, the height of the trend
#'   line does not seem robust.
#' 'MagTrendR': Relative difference between last value of trend prediction before
#'   and first value of trend prediction after drought breakpoint. Based on corrected
#'   trend for irregular data. Still, with irregular data, the height of the trend
#'   line does not seem robust.
#' 'AmpDiffR': Relative difference in mean amplitudes (based on sine and cosine terms
#'   of harmonic model) in segment before and after drought breakpoint.
#' 'Type': Typology of the drought breakpoint
#'
#' @details See \code{\link{resInd}}
#'  The function 'resIndSpatial' requires the function 'mc.calc' from the package
#'  'bfastSpatial' \url{https://github.com/loicdtx/bfastSpatial}, which you can install
#'  from github.
#'
#' @author Jennifer von Keyserlingk
#' @export
#'
#' @examples
#' #Load example raster data set (476 observations, 5x5 cells, including NA values)
#' #and date vector
#' data(stN) #raster brick
#' data(d) # date vector
#'
#' #Plot first 9 layers of raster brick with NDVI scaling factor
#' sc <- 0.0001
#' library(raster)
#' plot(stN*sc, 1:9)
#'
#' ##With package "bfastSpatial" you can extract information on acquisition date
#' # from typical Landsat file names \url{https://github.com/loicdtx/bfastSpatial}
#' #gs <- getSceneinfo(names(stN))
#' #d <- gs$date
#'
#'\dontrun{
#'#Run function "resIndSpatial". This requires the function 'mc.calc' from the
#' package "bfastSpatial", which you can install from github
#' # \url{https://github.com/loicdtx/bfastSpatial} .
#' y <- resIndSpatial(stN, d, dr=c(2004.753,2008.751), drd=2004.753)
#'
#'#Should return a raster stack containing 14 layers


resIndSpatial <- function(x, dates, type='irregular', sc=1, order=3,
                          formula = response ~ (trend + harmon), h=0.15,
                          plevel=0.05, dr, drd, s=3, NV=NA, mc.cores=1) {

fun <- function(x){

  ##Set output vector to NA
  resind <- NA

  ##Set own NoData value to differentiate between "no data pixels" and
  #"no breakpoint pixels"
  NV=NV

  #Apply bfastts() and multiply by data by scaling factor if required
  bfts <- bfastts(x*sc,dates,type=type)

  ##If not all oberservation NA (e.g. pixels outside area of interest) run procedure
  if(!all(is.na(bfts))){
    #Create data frame from bfastts
    bpp <- bfastpp(bfts, order=order, stl=("none"), na.action=na.omit)

    ##Calculate initial NDVI: First s years after start of observations
    Ini <- subset(bpp, time >= bpp$time[1] & time <= bpp$time[1]+s)
    MIni <- mean(Ini$response)

    ##MOSUM test for stuctural stability
    Mos <- efp(formula, data=bpp, type="OLS-MOSUM", h=h, rescale=FALSE)
    ##Calculate test statistics
    sct <- sctest(Mos)

    ##Fit breakpoints and segmented model if empirical fluctuation test was significant
    if (sct$p.value <= plevel){
      ##Fit the breakpoints
      #This is written in a tryCatch because sometimes an error in the RSS table occurs
      #It happens more in pixels with a smaller h and with many NA values
      #Probably due to singularities in the regression
      #However, it's unpredictable, so this tryCatch filters out those pixels
      #All procedures are terminated for these pixels and all indicators set to NA
      bpoints <- tryCatch({
        breakpoints(formula=formula, data=bpp, h=h)
      },
      error=function(cond) {
        message("Original error message:")
        message("Error in my.RSS.table[as.character(i), 3:4] <- c(pot.index[opt], break.RSS[opt]) :")
        message(cond[1])
        return(NA)
      }
      )
      # if the error did not occur and breakpoinits have been calculated, continue the procedure
      if(length(bpoints)==1){
        # if bpoints length=1 (NA) set bfts to NA too
        bfts <- NA
      } else {
        p <- bpoints$breakpoints
        if((length(p)==1 && is.na(p))) {
          #If no breakpoint p is NA. set breakpoint number to zero. Give warning.
          warning('No breakpoint found')
          bpnumb <- 0
          #Set breakpoint parameters to NA (needed for further calculations)
          bd <- NA
          cid <- NA
          #Fit unsegmented model
          m <- rlm (formula=formula, data=bpp, maxit=100)
          bpp$segment <- 1 ##even if you have only one segment you need to specify this,
          #for trend calculation later on
        } else {
          #If there is a breakpoint extract breakpoint parameters
          #breakpoint number
          bpnumb <- length(p)
          #breakdates
          bd <- bpp$time[p]
          #confidence invervals
          ci <- confint(bpoints, breakpoints=TRUE)
          bda <- ci[[1]]
          #cid <- bpp$time[c(bda[,1],bda[,3])]

          ##Fit segmented model
          #Create column "segment" that indices segment in the dataframe "bpp"
          bpp$segment <- breakfactor(bpoints)
          #RLM fit; I increased the default maxiterations (20 -> 100)
          m <- rlm (response ~ segment/(trend+harmon), data=bpp, maxit=100)
        }
      }
    } else {
      ##If MOSUM not significant set breakpoint number and DBP to zero and other
      #output variables to NV
      bpnumb <- 0

      #Fit unsegmented model
      m <- rlm (formula=formula, data=bpp, maxit=100)
      bpp$segment <- 1 ##even if you have only one segment you need to specify this,
      #for trend calculation later on
      warning('No significant deviation from structural stability (MOSUM
              test). No breakpoints fitted.')
    }
    if(!all(is.na(bfts))){

      #Predict values and add column "prediction" in in the dataframe "bpp"
      bpp$prediction <- predict(m,newdata=bpp)

      #Add trend prediction for each segment. Corrected hight of trend line for
      #irregular data: sets harmonic term based on mean DOY of observations/segment
      #(instead of trend$harmon[] <- 0, which assumes regular data);
      #idea based on email exchange with Achim Zeileis
      for(i in 1:length(unique(bpp$segment))) {
        seg <- unique(bpp$segment)[i]
        # trend <- subset(bpp, segment == seg)
        trend <- bpp[bpp$segment==seg, ]
        trend$days <- as.numeric(substr(formatC(trend$time,format='f',digits=3), 6, 8))
        dmean <- mean(trend$days)
        har <- dmean/365
        trend$harmon[] <- rep(
          c(cos(2 * pi * har * 1:order), sin(2 * pi * har * 1:order)),
          each = nrow(trend))

        #Add trendprediction column in bpp matrix
        bpp$trendprediction[bpp$segment == seg] <- predict(m, newdata = trend)
      }

      ##Add column for residuals to dataframe bpp
      # bpp$residuals <- as.numeric(bpp$response)-as.numeric(bpp$prediction)

      ##Extract coefficients
      coef <- coefficients(m)
      ##Save intercept
      Int <- coef[1]
      ##Extract trends
      trends <- coef[which(grepl('trend', names(coef)))]
      ##Extract Amplitudes for each segment. Depending on parameter "order" one
      # gets multiple sine and cosine terms
      cosines <- coef[which(grepl('harmoncos', names(coef)))]
      sines <- coef[which(grepl('harmonsin', names(coef)))]
      amps <- sqrt(cosines^2 + sines^2)
      names(amps) <- rep(1:length(unique(bpp$segment)), order)

      # Calculate (drought) resilience indicators -------------------------------

      ##If no (significant) breakpoint fuond set breakpoint position, breakpoint timing,
      #and recovery trend to "NV"; set DBP to 0
      if(bpnumb==0){
        DBP <- 0
        bpt <- NV
        tlag <- NV
        trendrecov <- NV
        pretrend <- NV
        preNDVI <- NV
        MagA <- NV
        MagR <- NV
        MagTA <- NV
        MagTR <- NV
        AmpDiff <- NV
        Type <- NV
        #Set breakpoint position to NA (needed for further calculations)
        bpd <- NA
      } else {
        ##If breakpoint, check if one breakpoint ocurred around drought period ("drought breakpoint")
        bpd <- match(bd[bd >= dr[1] & bd <= dr[2]], bd) #bpd gives you the position of bp.

        ##If drought breakpoint occured set DBP to 1 and calculate breakpoint timing,
        #trend of recovery, trend in segment before BP, preNDVI, Magnitude of change
        #(MagA & MagR) based on mean NDVI of s years before and after BP.
        #If there is more than 1 BP in the specified time interval, the first one is chosen!
        if (length(bpd) > 1) {
          bpd <- bpd[1]
        }

        if (length(bpd) > 0) {
          DBP <- 1

          #Calculate BP timing
          bpt <- bd[bpd]

          #Calculate tlag (days between drought reference year and bpt in days,
          # assuming 365 days/year)
          tlag <- (bpt-drd)*365

          #Calculate trend slopes
          seg2 <- (bpd+1)
          trendrecov <- trends[seg2]
          seg1 <- (bpd)
          pretrend <- trends[seg1]

          #Calculate preNDVI, Magnitude of change (MagA & MagR) based on mean NDVI
          #of s years before and after BP. Use subset of bpp$response that falls
          #within date vector
          ti <- c(bpt-s, bpt+s) # +/- s years around breakpoint
          S1 <- subset(bpp, time >= ti[1] & time <= bpt)
          S2 <- subset(bpp, time > bpt & time <= ti[2])

          preNDVI <- mean(S1$response)
          M2 <- mean(S2$response)

          MagA <- (M2-preNDVI) #absolute change magnitude
          MagR <- MagA/preNDVI #relative change magnitude

          #Calculate Magnitude of change based on corrected trend prediction (MagTA & MagTR)
          w <- which(bpp$time==bpt) #gives you position of breakpoint in trend
          #dataframe
          y1 <- bpp$trendprediction[w] #trend prediction value at time of breakpoint
          y2 <- bpp$trendprediction[w+1] #trend prediction value of observation after
          #breakpoint
          MagTA <- y2-y1 #absolute breakpoint magnitude
          MagTR <- (y2-y1)/y1 #relative breakpoint magnitude

          #Calculate Difference in mean amplitudes between segments
          mean_ampsseg1 <- mean(amps[which(grepl(seg1,names(amps)))]) #mean amplitude in segment before bp
          mean_ampsseg2 <- mean(amps[which(grepl(seg2,names(amps)))]) #mean amplitude in segment after bp

          AmpDiff <- (mean_ampsseg2-mean_ampsseg1)/mean_ampsseg1 #relative difference between mean amplitudes

          # Determine the breakpoint typology
          # interrupted increase (accelerating)
          if (pretrend>0 && trendrecov>0 && trendrecov>pretrend) {
            Type <- 1
          }
          # interrupted increase (slowing down)
          if (pretrend>0 && trendrecov>0 && trendrecov<pretrend) {
            Type <- 2
          }
          # interrupted decrease (accelerating)
          if (pretrend<0 && trendrecov<0 && trendrecov<pretrend) {
            Type <- 3
          }
          # interrupted decrease (slowing down)
          if (pretrend<0 && trendrecov<0 && trendrecov>pretrend) {
            Type <- 4
          }
          # positive reversal
          if (pretrend<0 && trendrecov>0) {
            Type <- 5
          }
          # negative reversal
          if (pretrend>0 && trendrecov<0) {
            Type <- 6
          }
        } else {
          ##If no breakpoint around drought occured set patrameters to NA
          warning('No drought breakpoint found')
          DBP <- 0
          bpt <- NV
          tlag <- NV
          trendrecov <- NV
          pretrend <- NV
          preNDVI <- NV
          MagA <- NV
          MagR <- NV
          MagTA <- NV
          MagTR <- NV
          AmpDiff <- NV
          Type <- NV
        }
      }
    }
    # If bpoints is NA due to error in bpoint calculation:
    if (length(bpoints)==1) {
      # If breakpoints was equal to NA
      warning('Error occurred during computation of breakpoints. All ouptut variables are set to NA')
      bpnumb <- NA
      MIni <- NA
      Int <- NA
      DBP <- NA
      bpt <- NA
      tlag <- NA
      trendrecov <- NA
      pretrend <- NA
      preNDVI <- NA
      MagA <- NA
      MagR <- NA
      MagTA <- NA
      MagTR <- NA
      AmpDiff <- NA
      Type <- NA
    }
  } else {
    #If all values NA (pixels without NDVI values!) set output variables to NA and give warning
    warning('No observations in time series')
    bpnumb <- NA
    MIni <- NA
    Int <- NA
    DBP <- NA
    bpt <- NA
    tlag <- NA
    trendrecov <- NA
    pretrend <- NA
    preNDVI <- NA
    MagA <- NA
    MagR <- NA
    MagTA <- NA
    MagTR <- NA
    AmpDiff <- NA
    Type <- NA
  }

  # Save and return output --------------------------------------------------

  #Save all indicators::
  resind <- cbind(bpnumb, MIni, as.numeric(Int), DBP, bpt, tlag, as.numeric(trendrecov),
                  as.numeric(pretrend), preNDVI, MagA, MagR, MagTA, MagTR, AmpDiff,Type)

  colnames(resind) <- c('BPNumb', 'Initial NDVI', 'Intercept', 'DBP','BpTime',
                        'Timelag', 'RecTrend', 'PreTrend', 'PreNDVI', 'MagObsA',
                        'MagObsR',   'MagTrendA', 'MagTrendR', 'AmpDiffR','Type')

  return(resind)
  }

  out <- bfastSpatial::mc.calc(x=x, fun=fun, mc.cores=mc.cores)
  names(out) <- c('BPNumb', 'Initial NDVI', 'Intercept', 'DBP','BpTime', 'Timelag',
                  'RecTrend', 'PreTrend', 'PreNDVI', 'MagObsA', 'MagObsR',
                  'MagTrendA', 'MagTrendR', 'AmpDiffR','Type')

  return(out)
}
