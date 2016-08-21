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
#'
#' @return
#'   list of reconciled time series. If \code{returnall==FALSE}, only the most 
#'   disaggregated series is returned.
#' @seealso \code{\link{thief}}, \code{\link{tsaggregates}}
#'
#' @examples
#' # Construct aggregates
#' aggts <- tsaggregates(USAccDeaths)
#' 
#' # Compute forecasts
#' fc <- list()
#' for(i in seq_along(aggts))
#'   fc[[i]] <- forecast(aggts[[i]], h=2*frequency(aggts[[i]]))$mean
#' names(fc) <- names(aggts)
#' 
#' # Reconcile forecasts
#' z <- reconcilethief(fc)
#' 
#' # Plot forecasts before and after reconcilation
#' par(mfrow=c(3,2))
#' for(i in seq_along(fc))
#' {
#'   plot(fc[[i]], ylab="", main=names(fc)[i],
#'     ylim=range(fc[[i]],z[[i]]))
#'   lines(z[[i]],col='red')
#' }
#'
#' @export
#' @author Rob J Hyndman

reconcilethief <- function(forecasts,
               comb=c("struc","mse","ols","bu","shr","sam"),
               mse=NULL, residuals=NULL, returnall=TRUE)
{
  comb <- match.arg(comb)

  # If forecasts is a list of forecast objects, then
  # extract list of forecast time series and list of residual time series
  if(is.element("forecast",class(forecasts[[1]])))
  {
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
    return(forecasts[[1]])

  tspf <- tsp(forecasts[[1]])

  # Set up group matrix for hts
  nsum <- rep(freq,m/freq)
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

  if(is.element(comb, c("struc","ols","mse"))) # OLS and WLS
  {
    if(comb=="struc")
      weights <- 1/nsum
    else if(comb=="ols")
      weights <- NULL
    else if(comb=="mse")
      weights <- 1/rep(rev(mse), rev(unsum))
    bts <- hts::combinef(t(fmat), groups=grps, weights=weights, keep='bottom')
  }
  else   # GLS
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
  bts <- ts(c(t(bts)))

  tsp(bts) <- tspf
  if(returnall)
    return(tsaggregates(bts))
  else
    return(bts)
}
