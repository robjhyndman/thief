#' Temporal hierarchical forecasting
#'
#' Takes a time series as input and produces forecasts using
#' the temporal hierarchical approach of Athanasopoulos et al (2016).
#'
#' @param y        Time series input
#' @param m        Seasonal period
#' @param h        Forecast horizon
#' @param comb     Combination method of temporal hierarchies, taking one of the following values:
#' \describe{
#'   \item{"struc"}{Structural scaling - weights from temporal hierarchy}
#'   \item{"mse"}{Variance scaling - weights from in-sample MSE}
#'   \item{"ols"}{Unscaled OLS combination weights}
#'   \item{"bu"}{Bottom-up combination -- i.e., all aggregate forecasts are ignored.}
#'   \item{"shr"}{GLS using a shrinkage (to block diagonal) estimate of residuals}
#'   \item{"sam"}{GLS using sample covariance matrix of residuals}
#' }
#' @param usemodel    Model used for forecasting each aggregation level:
#' \describe{
#'   \item{"ets"}{exponential smoothing, using the \code{\link[forecast]{ets}} function.}
#'   \item{"arima"}{arima, using the \code{\link[forecast]{auto.arima}} function.}
#'   \item{"theta"}{theta method, using the \code{\link[forecast]{thetaf}} function.}
#'   \item{"naive"}{random walk forecasts}
#'   \item{"snaive"}{seasonal naive forecasts, based on the last year of observed data.}
#' }
#' @param forecastfunction User-defined function to be used instead of \code{usemodel}. The
#' function must take a time series as the first argument, and the forecast horizon 
#' as the second argument. It must return an object of class \code{forecast}.
#' @param ...   Arguments to be passed to the time series modelling function 
#' (such as \code{ets} or \code{auto.arima}), or to \code{forecastfunction}.
#' 
#' @details This function computes the temporal aggregates of \code{y} using 
#' \code{\link{tsaggregates}}, then calculates all forecasts using the model function 
#' specified by \code{usemodel} or \code{forecastfunction}, and finally reconciles the 
#' forecasts using \code{\link{reconcilethief}}. The reconciled forecasts of \code{y} 
#' are returned.
#'
#' @return
#'   forecast object.
#' 
#' @seealso \code{\link{reconcilethief}}
#' 
#' @examples
#' \dontrun{
#' 
#' # Select ARIMA models for all series using auto.arima()
#' z <- thief(AEdemand[,12], usemodel='arima')
#' plot(z)
#' 
#' # Use your own function
#' ftbats <- function(y,h,...){forecast(tbats(y),h,...)}
#' z <- thief(AEdemand[,12], forecastfunction=ftbats)
#' plot(z)
#' }
#'
#' @export
#' @import forecast
#' @author Rob J Hyndman and Nikolaos Kourentzes

thief <- function(y, m=frequency(y), h=m*2,
               comb=c("struc","mse","ols","bu","shr","sam"),
               usemodel=c("ets","arima","theta","naive","snaive"), 
               forecastfunction=NULL, ...)
{
  comb <- match.arg(comb)
  if(is.null(forecastfunction))
    usemodel <- match.arg(usemodel)
  else
    usemodel <- deparse(substitute(forecastfunction))

  # Check input is a univariate time series
  if(!is.element("ts",class(y)))
    stop("y must be a time series object")
  if(NCOL(y) > 1L)
    stop("y must be a univariate time series")

  # Make sure the time series is seasonal
  if(m <= 1L)
    stop("Seasonal period (m) must be greater than 1")
  if(length(y) < 2*m)
    stop("I need at least 2 periods of data")

  # Compute aggregate
  aggy <- tsaggregates(y)

  # Compute forecasts
  frc <- th.forecast(aggy, h=h, usemodel=usemodel, 
    forecastfunction=forecastfunction, ...)

  # Reconcile forecasts and fitted values
  bts <- reconcilethief(frc$forecast, comb, frc$mse, frc$residuals, returnall=FALSE)
  fits <- reconcilethief(frc$fitted, comb, frc$mse, frc$residuals, returnall=FALSE)

  # Truncate to h forecasts
  tspb <- tsp(bts)
  bts <- ts(bts[1:h], start=tspb[1], frequency=tspb[3])

  # Construct forecast object to return
  out <- structure(list(x=y, mean=bts, fitted=fits, residuals=y-fits),
    class='forecast')
  out$method <- paste("THieF-",toupper(usemodel),sep="")
  return(out)
}

#-------------------------------------------------
th.forecast <- function(aggy, h=NULL, usemodel, forecastfunction, ...)
{
# Produce forecasts for multiple temporal aggregation levels with predefined models
#
# Inputs:
#   aggy        List of time series to be forecast
#   h           Forecast horizon. Default is m*2
#   usemodel    Model used for forecasting each aggregation level
#
# Outputs:
#   frc         Uncombined temporal hierarchy forecasts
#   fitted      Uncombined temporal hierarchy fitted values
#   residual    Uncombined temporal hierarchy residuals
#   mse         Mean squared errors for each temporal aggregation level
#   mseh        Mean squared errors per horizon for each temporal aggregation level
#
# Example:
#   aggy <- tsaggregates(USAccDeaths)
#   th.forecast(aggy,h=20)
#
# Nikolaos Kourentzes and Rob J Hyndman

  # Calculate forecast horizons
  n.k <- length(aggy)
  freq <- unlist(lapply(aggy, frequency))
  m <- max(freq)
  AL <- m/freq
  H <- m/AL*ceiling(h/m) 

  # Initialise
  frc <- fitted <- resid <- mseh <- vector("list",n.k)
  mse <- vector("numeric",n.k)

  # Model estimation and forecasts
  for (k in 1:n.k)
  {
    temp <- th.forecast.loop(aggy[[k]], H[[k]], usemodel, forecastfunction, ...)
    frc[[k]] <- temp$frc
    fitted[[k]] <- temp$fitted
    resid[[k]] <- temp$resid
    mse[k] <- temp$mse
    mseh[[k]] <- temp$mseh
  }

  # Return forecasts
  #names(frc) <- names(fitted) <- names(resid) <- names(mseh) <- names(mse) <- names(aggy)

  return(list(forecast=frc,fitted=fitted,residuals=resid,mse=mse,mseh=mseh))
}

#-------------------------------------------------
th.forecast.loop <- function(y, h, usemodel, forecastfunction, ...){
# Loop for model forecasts
# Produce forecasts and other information for univariate series y
# Forecast horizon h
# Return list of forecasts, fitted values, residuals, mse, etc.

  m <- frequency(y)

  if(!is.null(forecastfunction))
  {
    fc <- forecastfunction(y, h, ...)
    # Check if PI returned.
    if(!is.null(fc$lower))
    {
      # Check only 80% PI is saved
      if(is.element(80,fc$level))
        fc$lower <- fc$lower[,fc$level==80]
      else
        fc$lower <- NULL
    }
  }
  else if (usemodel=="ets")
  {
    fit <- try(forecast::ets(y, ...), silent=TRUE)
    # Do something simpler if ets model doesn't work
    # This is only necessary for versions of the package before 22 Aug 2016 (7.2beta)
    if(is.element("try-error",class(fit)))
    {
      # Use HW if there is enough data and the data is seasonal
      if(frequency(y) > 1L & length(y) >= 2*frequency(y))
        fit <- fc <- forecast::hw(y, h=h, level=80, 
                        initial='simple', alpha=0.2, beta=0.1, gamma=0.01)
      else # Otherwise just use Holt's
        fit <- fc <- forecast::holt(y, h=h, level=80, 
                        initial='simple', alpha=0.2, beta=0.1)
    }
    else
      fc <- forecast::forecast(fit, h=h, level=80)
  }
  else if(usemodel=="arima")
  {
    fit <- forecast::auto.arima(y, ...)
    fc <- forecast::forecast(fit, h=h, level=80)
  } 
  else if (usemodel == "theta")
  {
    fc <- forecast::thetaf(y, h=h, level=80)
  }
  else 
  {
    # Seasonal naive forecast
    if(usemodel=='snaive')
      fc <- forecast::snaive(y, h=h, level=80)
    # Naive forecast
    else
      fc <- forecast::naive(y, h=h, level=80)
  }

  # Calculate fitted residuals and MSE
  fitted <- fitted(fc)
  if(is.null(fitted) & !is.null(fc$residuals))
    fitted <- y - fc$residuals
  if(!is.null(fitted))
    resid <- y - fitted
  else
    resid <- NULL
  mse <- mean(resid^2, na.rm=TRUE)

  # Only keep full years of fitted values and residuals
  # Remove partial year from start
  tspy <- tsp(y)
  fullyears <- trunc(length(y) / m)
  fitted <- ts(utils::tail(fitted, fullyears*m), frequency=m, end=tspy[2])
  resid <- ts(utils::tail(resid, fullyears*m), frequency=m, end=tspy[2])

  # Compute MSE and MSEH
  if(!is.null(fc$lower))
    mseh <- ((fc$mean - fc$lower)/stats::qnorm(0.9))^2
  else
    mseh <- NULL

  return(list(frc=fc$mean,mse=mse,mseh=mseh,fitted=fitted,resid=resid))
}
