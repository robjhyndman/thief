#' Plot time series aggregates
#'
#' Plots all temporal aggregations of a time series
#'
#' @param x \code{tsaggregates} object, produced by \code{\link{tsaggregates}}.
#' @param object \code{tsaggregates} object, produced by \code{\link{tsaggregates}}.
#' @param series The indexes of the series to plot. By default, all series are plotted.
#' @param ... Other arguments passed to \code{\link[stats]{plot.ts}}
#' or \code{\link[forecast]{autoplot.ts}}.
#'
#' @method plot tsaggregates
#'
#' @examples
#' deathagg <- tsaggregates(USAccDeaths)
#' plot(deathagg, series=c(1,2,4,6))
#'
#' library(ggplot2)
#' autoplot(deathagg)
#'
#' @author Rob J Hyndman
#' @export
#' @importFrom stats plot.ts
#' @importFrom graphics plot

plot.tsaggregates <- function(x, series="all", ...)
{
  tmp <- mts.from.tsaggregates(x,series,...)
  args <- tmp$args
  args$x <- tmp$tsmat
  if(!is.element("main",names(args)))
    args$main <- deparse(substitute(x))
  if(length(tmp$series)==1L & !is.element("ylab",names(args)))
    args$ylab <- names(x)[tmp$series]
  if(length(tmp$series) <= 6L & !is.element("nc",names(args)))
    args$nc <- 1L
  do.call(plot.ts, args)
}

#' @importFrom ggplot2 autoplot
#' @rdname plot.tsaggregates
#' @method autoplot tsaggregates
#' @export

autoplot.tsaggregates <- function(object, series="all", ...)
{
  tmp <- mts.from.tsaggregates(object,series,...)
  args <- tmp$args
  args$object <- tmp$tsmat
  args$facets <- TRUE
  if(!is.element("main",names(args)))
    args$main <- deparse(substitute(object))
  if(!is.element("ylab",names(args)))
  {
    if(length(tmp$series)==1L)
      args$ylab <- names(object)[tmp$series]
    else
      args$ylab <- ""
  }
  do.call(ggplot2::autoplot, args)
}

#' @importFrom stats approx frequency is.ts time ts tsp "tsp<-"

mts.from.tsaggregates <- function(x, series, ...)
{
  nseries <- length(x)

  if(identical(series,"all"))
    series <- seq_len(nseries)
  if(length(series) > 10L)
    stop("Too many series to plot. Please select a subset.")
  if(min(series) < 1L | max(series) > nseries)
    stop(paste("series should be between 1 and",nseries))

  args <- list(...)

  # Create linearly interpolated mts object
  times <- time(x[[1]])
  tsmat <- matrix(0, nrow=length(times), ncol=nseries)

  tsmat[,1] <- x[[1]]
  if(nseries > 1L)
  {
    for(j in 2:nseries)
    {
      if(length(x[[j]]) > 1L)
        tsmat[,j] <- approx(time(x[[j]]), x[[j]], xout=times)$y
    }
  }
  tsmat <- ts(tsmat)
  tsp(tsmat) <- tsp(x[[1]])
  colnames(tsmat) <- names(x)
  tsmat <- tsmat[,series]

  if(!is.element("xlab",names(args)))
  {
    m <- frequency(x[[1]])
    if(m == 12L | m==4L | m==8760L | m==17520L)
      args$xlab <- "Year"
    else if(m==7L)
      args$xlab <- "Week"
    else if(m==24L | m==48L)
      args$xlab <- "Day"
    else if(m==168L | m==336L)
      args$xlab <- "Week"
    else
      args$xlab <- "Time"
  }

  return(list(tsmat=tsmat, args=args, series=series))
}
