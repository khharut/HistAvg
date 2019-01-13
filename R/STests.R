
#' Estimates Hurst exponent for a vector
#'
#' Calculates Hurst exponent, i.e. long-time memory for a vector
#' @author Harutyun Khachatryan
#' @param x numeric vector
#' @return Hurst exponent value mainly a number between zero and unity
#' @export
#' @importFrom stats lm
Hurst2 <- function(x)
{
  # half intervals of indices
  half <- function(N) sort(c(N, N[-length(N)] + ((diff(N) + 1) %/% 2)))
  # define the R/S scale
  rscalc <- function(x) {
    n <- length(x); y <- cumsum(x - mean(x))
    R <- diff(range(y)); S <- sd(x)
    return(R / S)
  }
  # set initial values
  X <- c(length(x))
  Y <- c(rscalc(x))
  N <- c(0, length(x) %/% 2, length(x))
  # compute averaged R/S for halved intervals
  while (min(diff(N)) >= 8) {
    xl <- c(); yl <- c()
    for (i in 2:length(N)) {
      rs <- rscalc(x[(N[i - 1] + 1):N[i]])
      xl <- c(xl, N[i] - N[i - 1])
      yl <- c(yl, rs)
    }
    X <- c(X, mean(xl))
    Y <- c(Y, mean(yl))
    # next step
    N <- half(N)
  }
  # apply linear regression
  rs_lm <- lm(log(Y) ~ log(X))
  return(unname(coefficients(rs_lm)[2]))
}

#' Estimates Hurst exponent for a vector, wrapper
#'
#' Wrapper for Hurst exponent calculator. Works on any type of input data.
#' @author Harutyun Khachatryan
#' @param x numeric vector
#' @return Hurst exponent value mainly a number between zero and unity
#' @export
HurstWrapper2 <- function(x)
{
  hvalue <- 0
  data <- na.omit(x)
  try(hexp <- Hurst2(data), silent = TRUE)
  if (!exists("hexp"))
  {
    hexp <- 0
  } else
  {
    if (is.na(hexp))
    {
      hexp <- 0
    }
  }
  hvalue <- hexp
  return(hvalue)
}

#' Estimates autoregression coefficient for a vector, wrapper
#'
#' Wrapper for autoregression coefficient calculator. Works on any type of input data.
#' @author Harutyun Khachatryan
#' @param x numeric vector
#' @return Autoregression coefficient value mainly a number between zero and unity
#' @export
#' @importFrom stats arima
ARIMAWrapper2 <- function(x)
{
  rho <- 0
  data <- na.omit(x)
  error.message <- try(arho <- unname(ar(data, order.max = 1, method = "mle")$ar[1]), silent = TRUE)
  if (!exists("arho"))
  {
    data <- scale(data, center = TRUE, scale = TRUE)
    error.message <- try(arho <- unname(arima(data, c(1, 0, 0), method = "CSS", optim.method = "BFGS")$coef[1]), silent = TRUE)
  } else
  {
    if (is.na(arho))
    {
      data <- scale(data, center = TRUE, scale = TRUE)
      error.message <- try(arho <- unname(arima(data, c(1, 0, 0), method = "CSS", optim.method = "BFGS")$coef[1]), silent = TRUE)
    }
  }
  if (!exists("arho"))
  {
    arho <- 0.92
  }
  else
  {
    if (is.na(arho))
    {
      arho <- 0.92
    }
  }
  rho <- arho
  return(rho)
}

#' Searches lost data around stationarity change points
#'
#' Checks whether stationarity change appears due to lost data
#' @author Harutyun Khachatryan
#' @param xtimes observation times of timeseries
#' @param obstimes observation times of raw data
#' @param lobs length of window around ind index
#' @return 0 when change point is not affected by lost data 0.1 otherwise
#' @export
LDTest2 <- function(xtimes, obstimes, lobs)
{
  idtime <- 24 * 60 * 60 / lobs
  cpvalue <- 0
  if (length(obstimes) > 1)
  {
    notime <- length(which((xtimes <= obstimes[length(obstimes)]) & (xtimes > obstimes[1])))
    ndtime <- (as.numeric(obstimes[length(obstimes)]) - as.numeric(obstimes[1])) / idtime
    if ((notime/ndtime) < 0.75) {cpvalue <- 0.1}
  }
  if (length(obstimes) == 1)
  {
    notime1 <- length(which((xtimes <= (obstimes[1] + 11 * idtime)) & (xtimes > obstimes[1])))
    notime2 <- length(which((xtimes <= obstimes[1]) & (xtimes > (obstimes[1] - 11 * idtime))))
    if ( (notime1 < 7) | (notime2 < 7) ) {cpvalue <- 0.1}
  }
  return(cpvalue)
}

#' Performs mean and variance change test
#'
#' Does mean and variance change test around index ind of time series with window length lobs
#' @author Harutyun Khachatryan, Tigran Khachikyan
#' @param itseries time series like vector
#' @param ind index around which test is performed
#' @param lobs length of window around ind index
#' @return 0 when nor mean neither variance change is detected,
#' 0.05 when mean or variance change is detected,
#' 0.1 when both mean and variance change is detected
#' @export
#' @importFrom stats qt
MVTest2 <- function(itseries, ind, lobs)
{
  cpvalue <- 0
  tstart <- max(c(ind - lobs + 1, 1))
  tend <- min(c(ind + lobs, length(itseries)))
  if ( (length(itseries[tstart:(ind - 1)]) > 1) & (length(itseries[(ind + 1):tend]) > 1) )
  {
    if ( (var(itseries[tstart:(ind - 1)], na.rm = TRUE) == 0) & (var(itseries[(ind + 1):tend], na.rm = TRUE) != 0))
    {
      cpvalue <- 0.1
    }
    if ( (var(itseries[tstart:(ind - 1)], na.rm = TRUE) != 0) & (var(itseries[(ind + 1):tend], na.rm = TRUE) == 0))
    {
      cpvalue <- 0.1
    }
    sindic <- abs(sd(itseries[tstart:(ind - 1)], na.rm = TRUE) - sd(itseries[(ind + 1):tend], na.rm = TRUE))
    sindic <- sindic * mean(itseries[tstart:tend], na.rm = TRUE)
    if (var(itseries[tstart:ind], na.rm = TRUE) * var(itseries[(ind + 1):tend], na.rm = TRUE) > 0)
    {
      intseries <- itseries[tstart:tend]
      rintseries <- rev(intseries)
      n <- length(intseries)
      mid_ind <- trunc(n / 2)
      v0 <- var(intseries[1:mid_ind], na.rm = TRUE)
      t_jn <- vector()
      rt_jn <- vector()
      v11 <- var(intseries[(n - mid_ind + 1):n], na.rm = TRUE)
      z <- v0 / v11
      rz <- v11 / v0
      if ((z > qf(1 - 0.001, mid_ind, mid_ind)) | (z < qf(0.001, mid_ind, mid_ind)) | (rz > qf(1 - 0.001, mid_ind, mid_ind)) | (rz < qf(0.001, mid_ind, mid_ind)))
      {
        cpvalue <- 0.05
      }
      for (j in 1:mid_ind)
      {
        v1 <- var(intseries[1:(j + mid_ind)], na.rm = TRUE)
        Sx_y <- sqrt( ((mid_ind * v0 + (j + mid_ind) * v1) *(2 * mid_ind + j)) / ((2 * mid_ind + j - 2) * mid_ind * (mid_ind + j)) )
        t_jn[j] <- (mean(intseries[1:mid_ind]) - mean(intseries[1:(j + mid_ind)])) / Sx_y
        rv1 <- var(rintseries[1:(j + mid_ind)], na.rm = TRUE)
        rSx_y <- sqrt(((mid_ind * v11 + (j + mid_ind) * rv1) * (2 * mid_ind + j)) / ((2 * mid_ind + j - 2) * mid_ind * (mid_ind + j)))
        rt_jn[j] <- (mean(rintseries[1:mid_ind]) - mean(rintseries[1:(j + mid_ind)])) / rSx_y
      }
      t_stat <- max(abs(c(t_jn, rt_jn)), na.rm = TRUE)
      dfree <- (n - 2)
      if (t_stat > qt(1 - 0.001, dfree, lower.tail = TRUE)) {cpvalue <- cpvalue + 0.05}
    }
    if (cpvalue > 0)
    {
      if (!is.na(sindic))
      {
        if (sindic <= 0.03) {cpvalue <- 0}
      }
    }
  }
  return(cpvalue)
}

#' Checks if any trend exists in data itseries
#'
#' Peroforms Spearman correlation test betwwen indeces and values of vector, if correlation exists then it must be some trend
#' @author Harutyun Khachatryan, Tigran Khachikyan
#' @param itseries time series like vector
#' @return 0 when nor mean neither variance change is detected,
#' 0.1 when both mean and variance change is detected
#' @export
STest2 <- function(itseries)
{
  cpvalue <- 0
  try(trho <- abs(cor.test(itseries, time(itseries), method = "spearman")$estimate), silent = TRUE)
  if (!exists("trho")) {trho <- 0} else {if (is.na(trho)) {trho <- 0}}
  if (trho > 0.9) {cpvalue <- 0.1}
  return(cpvalue)
}

#' Performs stationarity change test based on autoregression coefficient
#'
#' Does stationarity change test around index ind of time series with window length lobs
#' @author Harutyun Khachatryan, Tigran Khachikyan
#' @param itseries time series like vector
#' @param ind index around which test is performed
#' @param lobs length of window around ind index
#' @return 0 when nor mean neither variance change is detected,
#' 0.1 when both mean and variance change is detected
#' @export
ATest2 <- function(itseries, ind, lobs)
{
  cpvalue <- 0
  alpha <- 0.08
  arimastart <- max(c(ind - lobs + 1, 1))
  arimaend <- min(c(ind + lobs, length(itseries)))
  intseries <- as.vector(itseries[arimastart:arimaend])
  rho <- ARIMAWrapper2(intseries)
  if (rho >= (1 - alpha))
  {
    cpvalue <- 0.1
  }
  if (cpvalue > 0)
  {
    sindic <- abs(sd(itseries[seq(arimastart, ind - 1)], na.rm = TRUE) - sd(itseries[seq(ind + 1, arimaend)], na.rm = TRUE))
    sindic <- sindic * mean(itseries[arimastart:arimaend], na.rm = TRUE)
    if (!is.na(sindic))
    {
      if (sindic <= 0.03) {cpvalue <- 0}
    }
  }
  return(cpvalue)
}
