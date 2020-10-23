
#' @title Extracts change detection metrics from satellite time series
#' @description Computes several change detection metrics based on a BFAST (Break for Additive Season and Trend)
#' change detection framework. These change detection metrics are calculated for breakpoints occurring during a
#' given time period specified by the user (i.e. during a known drought). In addition, data about the full time
#' series is extracted: the overall number of breakpoints as well as the overall mean and initial value of the
#' dependent time series data. The included functions were designed for exploring the use of a breakpoint model
#' to study the response of vegetation to drought and climatic perturbartions. Please note that some of the
#' extracted change detection metrics are at an experimental stage. Applicable to individual pixels.
#'
#' @param x Numeric vector
#' @param dates Vector of dates in format "yyyy-mm-dd". Argument passed to bfastts()
#'  (see \code{\link[bfast]{bfastts}})
#' @param type Character. Time series type, "irregular", "16-day" or "10-day"
#'  (see \code{\link[bfast]{bfastts}}). Default="irregular"
#' @param sc Numeric. Scalefactor for input data. Default=1 (no rescaling).
#'  Scalefactor for input data. Default=1 (no rescaling).
#' @param order Numeric. Order of the harmonic component
#'  (see \code{\link[bfast]{bfastpp}}). Default=3.
#' @param formula Formula for the regression model and for fitting the breakpoints.
#'  Default=response ~ trend + harmon, a linear trend and a harmonic season component.#
#'  Argument passed to breakpoints() and efp() (see \code{\link[strucchange]{breakpoints}}
#'  and  \code{\link[strucchange]{efp}}).
#' @param h Minimal fraction of oberservations required between two breakpoints
#'  relative to the sample size. Default=0.15. Argument passed to efp() and to breakpoints()
#'  (see \code{\link[strucchange]{breakpoints}} and  \code{\link[strucchange]{efp}}).
#'  Default=0.15
#' @param plevel Numeric value. Critical significance level for MOSUM test of
#'  structural stability.
#'  If test result is above the critical value, no breakpoints will be fitted.
#'  Default=0.05
#' @param dr Date Vector c(min, max) setting the time period during which a breakpoint
#'  is searched. Format: c(decimal year, decimal year); assuming 365 days/year
#' @param drd Date. Reference day (e.g.start of drought) based on which the 'timelag'
#'  is calculated. Format: decimal years (assuming 365 days/year)
#' @param s Numeric. Number of years used to calculate mean NDVI after the beginning
#'  of the time series used for calculation mean 'initial NDVI', and number of years
#'  before and after the breakpoint used to calculate 'PreNDVI', 'MagObsA' and 'MagObsR'.
#'  Assumes 365 days/year. Default=3
#' @param NV Numeric. Sets the value assigned to pixels without a breakpoint during
#'  the specified time window. Default: NA
#' @param plot Logical. If TRUE a plot is produced. Default: FALSE
#'
#' @return numeric vector ('resind') with 14 entries:
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
#' @details This function was designed to explore several change detection metrics
#'  that can be extracted from a BFAST type change detection approach in relation
#'  to drought. Some of the extracted metrics are at the experimental stage and should be
#'  used with caution. Not all metrics are equally robust: for an irregular
#'  time series, the interecept of the linear trend line within segments (i.e. the
#'  height of the trend line after a breakpoint when plotted), was not stable.
#'  Therefore, the metrics relying directly on this information ('MagTrendA', 'MagTrendR')
#'  are equally not robust. However, robust results were obtainded for the slope of the
#'  trend line within segments ('RecTrend','PreTrend'), as well as the metrics 'BPBumb',
#'  'Initial NDVI', 'DBP', 'PreNDVI', 'MagObsR' and 'Type'.
#'  Those parts of the code dealing with the fitting of breakpoints to an irregular
#'  time series, as well as the fitting of BFAST type models to a segmented time series
#'  was based on the function "coefSegments" by Ben DeVries:
#'  \url{https://github.com/bendv/integrated-lts-cbm/blob/master/R/coefSegments.R}.
#'
#' @author Jennifer von Keyserlingk
#' @importFrom bfast bfastpp bfastts
#' @importFrom strucchange efp sctest breakpoints breakfactor
#' @importFrom MASS rlm
#' @importFrom zoo zoo
#' @importFrom graphics abline legend points lines
#' @importFrom stats coefficients confint na.omit time predict
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
#' #Select target pixel from raster stack
#' targcell <- 1
#'
#' #Extract vector of NDVI values for targcell and plot time series.
#' x <- as.vector(stN[targcell])
#' plot(d,x*sc, ylab='NDVI', xlab='year', main='targcell', ylim=c(0,1))
#'
#' #Run function "resInd"
#' y <- resInd(x, d, dr=c(2004.753,2008.751), drd=2004.753, plot=TRUE)
#'
#' #Should return a plot and a vector y containing 14 values


resInd <- function(x, dates, type='irregular', sc=1, order=3,
                   formula = response ~ (trend + harmon), h=0.15,
                   plevel=0.05, dr, drd, s=3, NV=NA, plot=FALSE) {
  #Set output vector to NA
  resind <- NA

  #Set own NoData value to differentiate between "no data pixels" and
  #"no breakpoint pixels"
  NV=NV

  #Apply bfastts() and multiply data by scaling factor if required
  bfts <- bfastts(x*sc,dates,type=type)

  ##If not all oberservation NA (e.g. pixels outside area of interest) run procedure
  if(!all(is.na(bfts))){
    #Create data frame from bfastts
    bpp <- bfastpp(bfts, order=order, stl=("none"), na.action=na.omit)

    #Calculate initial NDVI: First s years after start of observations
    Ini <- subset(bpp, time >= bpp$time[1] & time <= bpp$time[1]+s)
    MIni <- mean(Ini$response)

    ##MOSUM test for structural stability
    Mos <- efp(formula, data=bpp, type="OLS-MOSUM", h=h, rescale=FALSE)
    #Calculate test statistics
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
            #If no breakpoint p = NA. Set breakpoint number to zero. Give warning.
            warning('No breakpoint found')
            bpnumb <- 0
            #Set breakpoint parameters to NA (needed for further calculations)
            bd <- NA
            cid <- NA
            #Fit unsegmented model
            m <- rlm (formula=formula, data=bpp, maxit=100)
            bpp$segment <- 1 #even if you have only one segment you need to specify
            #this, for trend calculation later on
          } else {
            #If there is a breakpoint extract breakpoint parameters
            #breakpoint number
            bpnumb <- length(p)
            #breakdates
            bd <- bpp$time[p]
            #confidence invervals
            ci <- confint(bpoints, breakpoints=TRUE)
            bda <- ci[[1]]
            cid <- bpp$time[c(bda[,1],bda[,3])]

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
        warning('No significant deviation from structural stability (MOSUM test).
              No breakpoints fitted.')
    }
    if(!all(is.na(bfts))){

      #Predict values and add column "prediction" in in the dataframe "bpp"
      bpp$prediction <- predict(m,newdata=bpp)

      #Add trend prediction for each segment. Corrected hight of trend line for
      #irregular data: sets harmonic term based on mean DOY of observations/segment
      # (instead of trend$harmon[] <- 0, which assumes regular data);
      # idea based on email exchange with Achim Zeileis
      for(i in 1:length(unique(bpp$segment))) {
        seg <- unique(bpp$segment)[i]
        #trend <- subset(bpp, segment == seg)
        trend <- bpp[bpp$segment == seg,]
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
      ##Extract Amplitudes for each segment. Depending on parameter "order" one gets
      #multiple sine and cosine terms
      cosines <- coef[which(grepl('harmoncos', names(coef)))]
      sines <- coef[which(grepl('harmonsin', names(coef)))]
      amps <- sqrt(cosines^2 + sines^2)
      names(amps) <- rep(1:length(unique(bpp$segment)), order)

      # Calculate (drought) resilience indicators -------------------------------
      ##If no (significant) breakpoint was found set breakpoint position,
      # breakpoint timing, and recovery trend to "NV"; set DBP to 0
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
        ##If breakpoint, check if one breakpoint ocurred around drought period
        #("drought breakpoint")
        bpd <- match(bd[bd >= dr[1] & bd <= dr[2]], bd) #bpd gives you the position of bp.

        ##If drought breakpoint occured set DBP to 1 and calculate breakpoint timing,
        #trend of recovery, trend in segment before BP, preNDVI,
        #  Magnitude of change (MagA & MagR) based on mean NDVI of s years before and after BP.
        #If there is more than one BP in the specified time interval, the first one is chosen
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
          w <- which(bpp$time==bpt) #gives you position of breakpoint in trend dataframe
          y1 <- bpp$trendprediction[w] #trend prediction value at time of breakpoint
          y2 <- bpp$trendprediction[w+1] #trend prediction value of observation after breakpoint
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
    ##If all values NA (pixels without NDVI values!) set output variables to NA and give warning
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

  # Plotting ----------------------------------------------------------------

  if(!all(is.na(bfts))){
    if(plot) {

      #windows() to open plot in separate window

      plot(bpp$time,bpp$response, ylab="Data", xlab="Time", pch=4);
      points(bpp$time, bpp$prediction,col='red',type='b',pch=23);
      #Add vertical line at breakdates
      abline(v = bd, lty = 1, col="blue", lwd=1.2);
      #Add lines at dashed lines for confidence intervals
      abline(v = cid, lty = 3, col="blue", lwd=1.2);

      #Add legend
      legend("topleft",c("NDVI data", "fitted rlm season-trend model", "linear trend",
                         "breakdates", "confidence intervals breakdates"),
             col=c("black", "red", "grey28", "blue", "blue"), lty=c(0,1,1,1,3),
             pch=c(4,23,NA,NA,NA), cex=1, bty='n');

      #Add trend prediction for each segment. Corrected hight of trend line for
      #irregular data: fix harmonic based on mean DOY of observations/segment
      #(instead of trend$harmon[] <- 0, which assumes regular data);
      #idea based on email exchange with Achim Zeileis
      for(i in 1:length(unique(bpp$segment))) {
        seg <- unique(bpp$segment)[i]
        # trend <- subset(bpp, segment == seg)
        trend <- bpp[bpp$segment == seg,]
        trend$days <- as.numeric(substr(formatC(trend$time,format='f',digits=3), 6, 8))
        dmean <- mean(trend$days)
        har <- dmean/365
        trend$harmon[] <- rep(
          c(cos(2 * pi * har * 1:order), sin(2 * pi * har * 1:order)),
          each = nrow(trend))

        #Add trendprediction column in bpp matrix
        bpp$trendprediction[bpp$segment == seg] <- predict(m, newdata = trend)

        # plot trend
        lines(zoo::zoo(bpp$trendprediction[bpp$segment == seg],
                     bpp$time[bpp$segment == seg]), col = 'grey28')
      }
    } else {
      warning('No plot created')
      }
  }


  # Save and return output --------------------------------------------------

  #Save all indicators:
  resind <- cbind(bpnumb, MIni, as.numeric(Int), DBP, bpt, tlag,
                  as.numeric(trendrecov), as.numeric(pretrend), preNDVI, MagA,
                  MagR, MagTA, MagTR, AmpDiff,Type)

  colnames(resind) <- c('BPNumb', 'Initial NDVI', 'Intercept', 'DBP','BpTime',
                        'Timelag', 'RecTrend', 'PreTrend', 'PreNDVI', 'MagObsA',
                        'MagObsR', 'MagTrendA', 'MagTrendR', 'AmpDiffR','Type')

  return(resind)
}
