#' Reconcile temporal hierarchical forecasts
#'
#' Takes forecasts of time series at all levels of temporal aggregation
#' and combines them using the temporal hierarchical approach of Athanasopoulos et al (2016).
#'
#' @param forecasts List of forecasts. Each element must be a time series of forecasts,
#' or a forecast object.
#' The number of forecasts should be equal to k times the seasonal period for each series,
#' where k is the same across all series.
#' @param comb Combination method of temporal hierarchies, taking one of the following values:
#' \describe{
#'   \item{"struc"}{Structural scaling - weights from temporal hierarchy}
#'   \item{"mse"}{Variance scaling - weights from in-sample MSE}
#'   \item{"ols"}{Unscaled OLS combination weights}
#'   \item{"bu"}{Bottom-up combination -- i.e., all aggregate forecasts are ignored.}
#'   \item{"shr"}{GLS using a shrinkage (to block diagonal) estimate of residuals}
#'   \item{"sam"}{GLS using sample covariance matrix of residuals}
#' }
#' @param mse A vector of one-step MSE values corresponding to each of the forecast series.
#' @param residuals List of residuals corresponding to each of the forecast models.
#' Each element must be a time series of residuals. If \code{forecast} contains a list of
#' forecast objects, then the residuals will be extracted automatically and this argument
#' is not needed. However, it will be used if not \code{NULL}.
#' @param returnall If \code{TRUE}, a list of time series corresponding to the first argument
#' is returned, but now reconciled. Otherwise, only the most disaggregated series is returned.
#' @param aggregatelist (optional) User-selected list of forecast aggregates to consider
#'
#' @return
#'   List of reconciled forecasts in the same format as \code{forecast}.
#' If \code{returnall==FALSE}, only the most disaggregated series is returned.
#' @seealso \code{\link{thief}}, \code{\link{tsaggregates}}
#'
#' @examples
#' # Construct aggregates
#' aggts <- tsaggregates(USAccDeaths)
#'
#' # Compute forecasts
#' fc <- list()
#' for(i in seq_along(aggts))
#'   fc[[i]] <- forecast(auto.arima(aggts[[i]]), h=2*frequency(aggts[[i]]))
#'
#' # Reconcile forecasts
#' reconciled <- reconcilethief(fc)
#'
#' # Plot forecasts before and after reconcilation
#' par(mfrow=c(2,3))
#' for(i in seq_along(fc))
#' {
#'   plot(reconciled[[i]], main=names(aggts)[i])
#'   lines(fc[[i]]$mean, col='red')
#' }
#'
#' @export
#' @author Rob J Hyndman

reconcilethief <- function(forecasts,
               comb=c("struc","mse","ols","bu","shr","sam"),
               mse=NULL, residuals=NULL, returnall=TRUE, 
               aggregatelist=NULL)
{
  comb <- match.arg(comb)

  # If forecasts is a list of forecast objects, then
  # extract list of forecast time series and list of residual time series
  if(is.element("forecast",class(forecasts[[1]])))
  {
    returnclass <- "forecast"
    origf <- forecasts
    # Grab residuals
    if(is.null(residuals))
    {
      residuals <- list()
      for(i in seq_along(forecasts))
        residuals[[i]] <- residuals(forecasts[[i]])
      # Discard partial years at start of residual series
      for(i in seq_along(residuals))
      {
        tspy <- tsp(residuals[[i]])
        m <- tspy[3]
        fullyears <- trunc(length(residuals[[i]])/m)
        residuals[[i]] <- ts(utils::tail(residuals[[i]], fullyears*m), frequency=m, end=tspy[2])
      }
    }
    # Grab forecasts
    for(i in seq_along(forecasts))
      forecasts[[i]] <- forecasts[[i]]$mean
  }
  else
    returnclass <- "ts"

  # Find seasonal periods
  freq <- unlist(lapply(forecasts,frequency))
  if(min(freq) != 1L)
    stop("Minimum seasonal period should be 1")
  m <- max(freq)
  # Put in order
  k <- rev(order(freq))
  forecasts <- forecasts[k]
  freq <- freq[k]

  # Check series each have same equivalent lengths.
  lengths <- unlist(lapply(forecasts, length))
  if(!is.constant(lengths/freq))
    stop("Forecast series must have the same equivalent lengths")

  # Bottom up
  if(comb == "bu")
    bts <- forecasts[[1]]
  # Some combination
  else
  {
    # Set up group matrix for hts
    # (adjusted to allow consideration of aggregatelist input)
    nsum <- rev(rep(m/freq, freq))
    unsum <- unique(nsum)
    grps <- matrix(0, nrow=length(unsum)-1, ncol=m)
    for(i in 1:(length(unsum)-1))
    {
      mi <- m/unsum[i]
      grps[i,] <- rep(1:mi, rep(unsum[i],mi))
    }

    # Set up matrix of forecasts in right structure
    nc <- length(forecasts[[1]])/m
    fmat <- matrix(0, nrow=0, ncol=nc)
    for(i in rev(seq_along(forecasts)))
      fmat <- rbind(fmat, matrix(forecasts[[i]], ncol=nc))

    # OLS or WLS reconciliation
    if(is.element(comb, c("struc","ols","mse")))
    {
      if(comb=="struc")
        weights <- 1/nsum
      else if(comb=="ols")
        weights <- NULL
      else if(comb=="mse")
        weights <- 1/rep(rev(mse), rev(unsum))
      bts <- hts::combinef(t(fmat), groups=grps, weights=weights, keep='bottom')
    }
    else   # GLS reconciliation
    {
      # Set up matrix of residuals in right structure
      if(is.null(residuals))
        stop("GLS needs residuals")
      nc <- length(residuals[[1]])/m
      rmat <- matrix(0, nrow=0, ncol=nc)
      for(i in rev(seq_along(forecasts)))
        rmat <- rbind(rmat, matrix(residuals[[i]], ncol=nc))
      bts <- hts::MinT(t(fmat), groups=grps, residual=t(rmat),  covariance=comb, keep='bottom')
    }
    # Turn resulting reconciled forecasts back into a ts object
    bts <- ts(c(t(bts)))
    tsp(bts) <- tsp(forecasts[[1]])
  }

  # Now figure out what to return
  if(returnclass=="ts") # Just return time series
  {
    if(!returnall)
      return(bts)
    else
      return(tsaggregates(bts, aggregatelist=aggregatelist))
  }
  else #return forecast objects
  {
    if(!returnall)
    {
      adj <- bts - origf[[1]]$mean
      origf[[1]]$mean <- bts
      if(!is.null(origf[[1]]$lower))
      {
        origf[[1]]$lower <- sweep(origf[[1]]$lower,1,adj,"+")
        origf[[1]]$upper <- sweep(origf[[1]]$upper,1,adj,"+")
      }
      return(origf[[1]])
    }
    else
    {
      allts <- tsaggregates(bts, aggregatelist=aggregatelist)
      for(i in seq_along(origf))
      {
        adj <- allts[[i]] - origf[[i]]$mean
        origf[[i]]$mean <- allts[[i]]
        if(!is.null(origf[[i]]$lower))
        {
          origf[[i]]$lower <- sweep(origf[[i]]$lower,1,adj,"+")
          origf[[i]]$upper <- sweep(origf[[i]]$upper,1,adj,"+")
        }
      }
      return(origf)
    }
  }
}
