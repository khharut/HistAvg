
#' Calculates main periodicity markers for itseries
#'
#' Estimates main statistical parameters that shows that data is periodic
#' @author Harutyun Khachatryan, Tigran Khachikyan
#' @param itseries timeseries like vector of values
#' @param nfreq number of steps of possible periodicity to be checked
#' @return vector of 4 statistical parameters showing periodicity,
#' a) last period seasonal component and remainder component standard deviation ratio,
#' b) Spearman correlation rho for last period and one before last one
#' c) Spearman correlation rho for last period and mean periodic pattern of data
#' d) median value of seasonal and remainder components ratio for last period
#' @export
DRel2 <- function(itseries, nfreq)
{
  ratio <- 0
  relat <- 0
  crho_last <- 0
  coeff <- 1
  intseries <- na.omit(itseries)
  intseries <- intseries - min(intseries, na.rm = TRUE) + 1
  nobsi <- length(intseries)
  crho <- 0
  if ((nobsi / nfreq) > 2)
  {
    indi <- seq((nobsi - nfreq + 1), nobsi)
    intseries_right <- intseries[indi]
    intseries_left <- intseries[(indi - nfreq)]
    error.message <- try(crho <- suppressWarnings(cor.test(intseries_left, intseries_right, method = "spearman", exact = FALSE)$estimate), silent = TRUE)
    if (!exists("crho")) {crho <- 0} else {crho <- unname(crho); if (is.na(crho)) {crho <- 0}}
  }
  if ((nobsi / nfreq) >= 3)
  {
    xad.ts <- ts(intseries, frequency = nfreq)
    error.message <- try(xad.decomp <- stl(xad.ts, s.window = "per", robust = FALSE), silent = TRUE)
    if (exists("xad.decomp"))
    {
      ses <- xad.decomp$time.series[indi, "seasonal"]
      rem <- xad.decomp$time.series[indi, "remainder"]
      if (sd(rem, na.rm = TRUE) > 0)
      {
        error.message <- try(crho_last <- suppressWarnings(cor.test(ses, intseries_right, method = "spearman", exact = FALSE)$estimate), silent = TRUE)
        if (!exists("crho_last")) {crho_last <- 0} else {crho_last <- unname(crho_last); if (is.na(crho_last)) {crho_last <- 0}}
        ratio <- sd(ses, na.rm = TRUE) / sd(rem, na.rm = TRUE)
        relat <- median(ses / rem, na.rm = TRUE)
        if ( (relat < 0)  & (crho < 0.4) & (crho_last < 0.4) ) {coeff <- 0}
        if ( (min(c(crho,crho_last)) < 0.1) & (relat < 0) )  {coeff <- 0}
        if ( (crho > 0.7) | (crho_last > 0.9) ) {coeff <- 2}
      } else
      {
        if (sd(ses, na.rm = TRUE) > 0)
        {
        	ratio <- Inf
        	relat <- Inf
        }
      }
    }
  }
  ratio <- ratio * coeff
  return(c(ratio, crho, crho_last, relat, coeff))
}

#' Checks periodicity and gives reliability for having periodicity
#'
#' Performs different tests including ACF to detect periodicity in data
#' @author Harutyun Khachatryan, Tigran Khachikyan
#' @param itseries timeseries like vector of values
#' @param obsi number of steps of possible periodicity to be checked
#' @return vector of 2 elements containing:
#' a) input value obsi
#' b) reliability value spans from zero to unity
#' @export
#' @importFrom stats acf
#' @importFrom simecol peaks
NRel2 <- function(itseries, obsi)
{
  nacfpeak <- c(obsi, 0)
  intseries <- as.vector(na.omit(itseries))
  instep <- length(intseries)
  if ((instep / obsi) > 2)
  {
    x1d.acf <- acf(intseries, lag.max = ceiling(1.5 * obsi), plot = FALSE)
    x1d.acf.lag <- x1d.acf$lag
    x1d.acf.acf <- x1d.acf$acf
    x1d.acf.lag <- x1d.acf.lag[-1]
    x1d.acf.acf <- x1d.acf.acf[-1]
    peak1d <- peaks(x1d.acf.lag, x1d.acf.acf, mode = "max")
    if (length(which( (peak1d$x >= (obsi - (obsi / 12))) & (peak1d$x <= (obsi + (obsi/12))) ) ) > 0)
    {
      corval <- max(peak1d$y[which((peak1d$x >= (obsi - (obsi / 12))) & (peak1d$x <= (obsi + (obsi / 12))))], na.rm = TRUE)
    } else {corval <- 0}
    per_markers <- DRel2(intseries, obsi)
    if ((corval) < 0.1)
    {
      tin <- 0
    } else
    {
      tin <- 1
    }
    if (corval < 0.25)
    {
      corval <- 0
      if ((per_markers[2] < 0.1) | (per_markers[3] < 0.5)) {tin <- 0}
    }
    reliab = corval + tin * per_markers[1]
    reliab = atan(reliab) / (2 * atan(1))
    if ( ((instep / obsi) < 3) & (per_markers[1] == 0) )
    {
      if ( (corval > 0.25) & (per_markers[2] > 0.5) ) {reliab = max(c(corval, per_markers[2]))}
    }
    nacfpeak <- c(obsi, reliab)
  }
  return(nacfpeak)
}
