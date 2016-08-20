#' Temporal hierarchy forecasting
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
#' @param model    Model used for forecasting each aggregation level:
#' \describe{
#'   \item{"ets"}{exponential smoothing, using the \code{\link[forecast]{ets}} function.}
#'   \item{"arima"}{arima, using the \code{\link[forecast]{auto.arima}} function.}
#'   \item{"theta"}{theta method, using the \code{\link[TStools]{theta}} function.}
#'   \item{"naive"}{random walk forecasts}
#'   \item{"snaive"}{seasonal naive forecasts, based on the last year of observed data.}
#' }
#' @param forecastfunction User-defined function to be used instead of \code{model}. The
#' function must take a time series as the first argument, and the forecast horizon 
#' as the second argument. It must return an object of class \code{forecast}.
#' @param ...   Arguments to be passed to the time series modelling function 
#' (such as \code{ets} or \code{auto.arima}), or to \code{forecastfunction}.
#'
#' @return
#'   forecast object for bottom level series.
#'
#' @examples
#'   z <- thief(AEdemand[,12], model='arima')
#'   plot(z)
#'
#' @export
#' @import forecast
#' @author Rob J Hyndman and Nikolaos Kourentzes

thief <- function(y, m=frequency(y), h=m*2,
               comb=c("struc","mse","ols","bu","shr","sam"),
               model=c("ets","arima","theta","naive","snaive"), 
               forecastfunction=NULL, ...)
{
  comb <- match.arg(comb)
  model <- match.arg(model)

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
  aggy <- tsaggregates(y=y, m=m, align='end')

  # Compute forecasts
  frc <- th.forecast(aggy, h=h, m=m, model=model, 
    forecastfunction=forecastfunction, ...)

  # Set up group matrix for hts
  nsum <- as.numeric(sub("AL","",rownames(frc$forecast)))
  unsum <- unique(nsum)
  grps <- matrix(0, nrow=length(unsum)-1, ncol=m)
  for(i in 1:(length(unsum)-1))
  {
    mi <- m/unsum[i]
    grps[i,] <- rep(1:mi, rep(unsum[i],mi))
  }

  # Bottom up
  if(comb == "bu")
  {
    bts <- ts(c(utils::tail(frc$forecast, m)), frequency=m, start=tsp(y)[2] + 1/m)
    fits <- ts(c(utils::tail(frc$fitted, m)), frequency=m, start=tsp(y)[1])
  }
  else 
  {
    # OLS and WLS
    if(is.element(comb, c("struc","ols","mse")))
    {
      if(comb=="struc")
        weights <- 1/nsum
      else if(comb=="ols")
        weights <- NULL
      else if(comb=="mse")
        weights <- 1/rep(rev(frc$mse), rev(unsum))
      bts <- hts::combinef(t(frc$forecast), groups=grps, weights=weights, keep='bottom')
      fits <- hts::combinef(t(frc$fitted), groups=grps, weights=weights, keep='bottom')
    }
    # GLS
    else
    {
      bts <- hts::MinT(t(frc$forecast), groups=grps, residual=t(frc$residual), 
        covariance=comb, keep='bottom')
      fits <- hts::MinT(t(frc$fitted), groups=grps, residual=t(frc$residual), 
        covariance=comb, keep='bottom')
    }
    # Save bottom level forecasts
    bts <- ts(c(t(bts)), frequency=m, start=tsp(y)[2] + 1/m)
    # Reconcile fitted values
    fits <- ts(c(t(fits)), frequency=m, start=tsp(y)[1])
  }
  
  out <- structure(list(x=y, mean=bts, fitted=fits, residuals=y-fits),
    class='forecast')
  out$method <- paste("THieF-",toupper(model),sep="")
  return(out)
}

#-------------------------------------------------
th.forecast <- function(aggy, h=NULL, m, model, forecastfunction, thr=3, ...)
{
# Produce forecasts for multiple temporal aggregation levels with predefined models
#
# Inputs:
#   aggy        List of time series to be forecast
#   h           Forecast horizon. Default is m*2
#   model       Model used for forecasting each aggregation level
#   thr         Threshold for identifying outliers for theta.out. Outliers are standard errors >= thr.
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

  # Produce forecasts
  # Calculate maximum horizon
  n.k <- length(aggy)
  AL <- m/unlist(lapply(aggy, frequency))
  H <- m/AL*ceiling(h/m)  # Forecast horizon
  # Initialise
  frc <- fitted <- resid <- mseh <- vector("list",n.k)
  mse <- vector("numeric",n.k)

  # Model estimation and forecasts
  for (k in 1:n.k)
  {
    temp <- th.forecast.loop(k,aggy,m,AL,H,model,thr, forecastfunction, ...)
    frc[[k]] <- temp$frc
    fitted[[k]] <- temp$fitted
    resid[[k]] <- temp$resid
    mse[k] <- temp$mse
    mseh[[k]] <- temp$mseh
  }

  names(frc) <- paste0("AL",AL)
  names(mse) <- paste0("AL",AL)

  # Re-order aggregation levels
  frc <- frc[n.k:1]
  # Create forecast vector across aggregation levels
  frc <- do.call("rbind",frc)
  # Reformat it as matrix in case it was a single column before
  frc <- matrix(frc,ncol=ceiling(h/m))

  rownames(frc) <- paste0("AL",rep(AL[n.k:1],m/AL[n.k:1]))
  colnames(frc) <- paste0("AL",m,"(t+",1:ceiling(h/m),")")

  fitted <- fitted[n.k:1]
  resid <- resid[n.k:1]
  n.in <- unlist(lapply(fitted,length))
  temp.fitted <- temp.resid <- array(NA,c(sum(m/AL),n.in[1]/(m/AL[n.k])))
  # Start/stop of each aggregation level
  idx <- c(1,cumsum((m/AL)[n.k:1])+1)
  for (i in 1:(n.in[1]/(m/AL[n.k]))) # For each full year
  {
    for (k in 1:n.k) # For each aggregation level
    {
      # Start from end towards start
      temp.fitted[idx[k]:(idx[k+1]-1), ((n.in[1]/(m/AL[n.k]))-i+1)] <-
        fitted[[k]][(n.in[k]-(m/AL[n.k-k+1])*i+1):(n.in[k]-(m/AL[n.k-k+1])*(i-1))]
      temp.resid[idx[k]:(idx[k+1]-1), ((n.in[1]/(m/AL[n.k]))-i+1)] <-
        resid[[k]][(n.in[k]-(m/AL[n.k-k+1])*i+1):(n.in[k]-(m/AL[n.k-k+1])*(i-1))]
    }
  }
  fitted <- matrix(temp.fitted,nrow=sum(m/AL))
  resid <- matrix(temp.resid,nrow=sum(m/AL))
  rownames(fitted) <- rownames(resid) <- rownames(frc)
  
  # MSEH
  if(length(mseh) == n.k)
  {
    names(mseh) <- paste0("AL",AL)
    mseh <- mseh[n.k:1]
    # Create forecast vector across aggregation levels
    mseh <- do.call("rbind",mseh)
    # Reformat it as matrix in case it was a single column before
    mseh <- matrix(mseh,ncol=ceiling(h/m))
    rownames(mseh) <- paste0("AL",rep(AL[n.k:1],m/AL[n.k:1]))
    colnames(mseh) <- paste0("AL",m,"(t+",1:ceiling(h/m),")")
  }
  else
    mseh <- NULL

  # Return forecasts
  return(list(forecast=frc,fitted=fitted,residual=resid,mse=mse,mseh=mseh))
}

#-------------------------------------------------
th.forecast.loop <- function(k,Y,m,AL,H,model,thr, forecastfunction, ...){
# Parallel loop for model forecasts
# Internal function

  mseh <- fitted <- resid <- NULL

  #print(paste("Fitting series of frequency",frequency(Y[[k]])))

  if(!is.null(forecastfunction))
  {
    fc <- forecastfunction(Y[[k]], H[k], ...)
    frc <- matrix(fc$mean, nrow=m/AL[k])
    fitted <- fitted(fc)
    if(is.null(fitted) & !is.null(fc$residuals))
      fitted <- Y[[k]] - fc$residuals
    if(!is.null(fc$lower))
    {
      mseh <- matrix(((fc$mean - fc$lower[,1])/abs(stats::qnorm((1-0.8)/2)))^2,
          nrow=m/AL[k])
    }
  }
  else if (model=="ets")
  {
    fit <- try(forecast::ets(Y[[k]], ...), silent=TRUE)
    # Do something simpler if ets model doesn't work
    if(is.element("try-error",class(fit)))
    {
      # Use HW if there is enough data and the data is seasonal
      if(frequency(Y[[k]]) > 1L & length(Y[[k]]) >= 2*frequency(Y[[k]]))
      {
        fit <- temp.frc <- forecast::hw(Y[[k]], h=H[k], level=95, 
                            initial='simple', alpha=0.2, beta=0.1, gamma=0.01)
      }
      else # Otherwise just use Holt's
      {
        fit <- temp.frc <- forecast::holt(Y[[k]], h=H[k], level=95, 
                               initial='simple', alpha=0.2, beta=0.1)
      }
    }
    else
    {
      temp.frc <- forecast::forecast(fit, h=H[k], level=95)
    }
    fitted <- fit$fitted
    frc <- matrix(temp.frc$mean,nrow=m/AL[k])
    mseh <- matrix(((temp.frc$mean - temp.frc$lower[,1])/abs(stats::qnorm((1-0.8)/2)))^2,
        nrow=m/AL[k])
  }
  else if(model=="arima")
  {
    fit <- forecast::auto.arima(Y[[k]], ...)
    temp.frc <- forecast::forecast(fit, h=H[k], level=95)
    fitted <- fit$x - fit$residuals
    frc <- matrix(temp.frc$mean, nrow=m/AL[k])
    mseh <- matrix(((temp.frc$mean - temp.frc$lower[,1])/abs(stats::qnorm((1-0.8)/2)))^2,
        nrow=m/AL[k])
  } 
  else if (model == "theta" || model == "theta.out")
  {
    # Theta
    cma <- TStools::cmav(Y[[k]])
    loc <- NULL
    if ((m/AL[k]) != 1)
    {
      is.mult <- TStools::mseastest(Y[[k]],cma=cma)$is.multiplicative
      # Identify outliers using time series decomposition
      if (model == "theta.out")
      {
        y.dcp <- TStools::decomp(Y[[k]],trend=cma)
        loc <- TStools::residout(y.dcp$irregular,t=thr)$location
      }
    } 
    else 
    {
      is.mult <- TRUE
      # Identify outliers using time series decomposition
      if (model == "theta.out")
        loc <- TStools::residout(Y[[k]]-cma,t=thr)$location
    }

    fit <- TStools::theta(Y[[k]],h=H[k], multiplicative=is.mult, outliers=loc)
    frc <- matrix(as.numeric(fit$frc),nrow=m/AL[k])
    fitted <- fit$fit
  }
  else 
  {
    # Seasonal naive forecast
    if(model=='snaive')
      fit <- forecast::snaive(Y[[k]], h=H[k], level=95)
    # Naive forecast
    else
      fit <- forecast::naive(Y[[k]], h=H[k], level=95)
    fitted <- fit$x - fit$residuals
    frc <- matrix(fit$mean, nrow=m/AL[k])
    mseh <- matrix(((fit$mean - fit$lower[,1])/abs(stats::qnorm((1-0.8)/2)))^2,
        nrow=m/AL[k])
  }

  # Calculate fitted residuals and MSE
  resid <- Y[[k]] - fitted
  mse <- mean(resid^2, na.rm=TRUE)

  return(list(frc=frc,mse=mse,mseh=mseh,fitted=fitted,resid=resid))
}
