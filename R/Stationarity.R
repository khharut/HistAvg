#' Single stationarity change
#'
#' Estimates single most possible stationarity point location in data
#' @author Harutyun Khachatryan
#' @param data input data vector
#' @return most possible change point location index
#' @export
#' @examples x<- rnorm(1000)
#' SingleChange2(x)
SingleChange2 <- function(data)
{
  n <- length(data)
  y <- diffinv(data)
  y2 <- diffinv((data) ^ 2)
  taustar <- 2:(n - 2)
  sigma1 <- ((y2[taustar + 1] - (y[taustar + 1] ^ 2/taustar))/taustar)
  neg <- sigma1 <= 0
  sigma1[neg == TRUE] <- 1 * 10 ^ (-10)
  sigman <-  ((y2[n + 1] - y2[taustar + 1]) - ((y[n + 1] -  y[taustar + 1]) ^ 2/(n - taustar)))/(n - taustar)
  neg <- sigman <= 0
  sigman[neg == TRUE] <- 1 * 10 ^ (-10)
  tmp <- taustar * log(sigma1) + (n - taustar) * log(sigman)
  tau <- which(tmp == min(tmp, na.rm = T))[1]
  tau <- tau + 1
  return(tau)
}

#' Multiple stationarity change
#'
#' Estimates possible stationarity point location in data which number at least is equal to Q
#' @author Harutyun Khachatryan
#' @param intseries input data vector
#' @param Q maximal number of changepoints to be found
#' @param robust logical TRUE/FALSE indicating if robust method should be used
#' @return most possible change point location index
#' @export
#' @examples x<- rnorm(1000)
#' MultipleChange2(x, Q = 5)
MultipleChange2 <- function(intseries, Q = 5, robust = FALSE)
{
  n <- length(intseries)
  if (n < 4)
  {
    stop("Data must have atleast 4 observations to fit a changepoint model.")
  }
  if (Q > ((n / 2) + 1))
  {
    stop(paste("Q is larger than the maximum number of segments", (n / 2) + 1))
  }
  std_data <- sd(intseries, na.rm = TRUE)
  mean_data <- mean(intseries, na.rm = TRUE)
  cpts <- NULL
  if (is.na(std_data))
  {
    std_data <- 0
  } else
  {
    if ((!is.numeric(std_data)) | (!is.numeric(mean_data)))
    {
      std_data <- 0
    }
  }
  if (std_data > 0)
  {
    if (robust)
    {
      data <- intseries - mean_data
    } else
    {
      data <- intseries / std_data
    }
    pen <- log(log(n))
    y2 <- c(0, cumsum(data ^ 2))
    tau <- c(0, n)
    cpt <- matrix(0, nrow = 2, ncol = Q)
    oldmax <- Inf
    lambda <- rep(0, n - 1)
    for (q in 1:Q)
    {
      i <- 1
      st <- tau[1] + 1
      end <-  tau[2]
      for (j in 1:(n - 1))
      {
        if (j == end)
        {
          st <- end + 1
          i <- i + 1
          end <- tau[i + 1]
        }
        else
        {
          if (y2[end + 1] != y2[st])
          {
            lambda[j] <- sqrt((end - st + 1)/2) * ((y2[j + 1] - y2[st]) / (y2[end + 1] - y2[st]) - (j - st + 1) / (end - st + 1))
          }
        }
      }
      k = which.max(abs(lambda))
      cpt[1, q] <- k
      cpt[2, q] <- min(oldmax, max(abs(lambda), na.rm = TRUE))
      oldmax <- min(oldmax, max(abs(lambda), na.rm = TRUE))
      tau <- sort(c(tau, k))
      lambda[] = 0
    }
    op.cps <- NULL
    criterion <- (cpt[2, ]) >= pen
    if (sum(criterion) == 0)
    {
      op.cps <- 0
    }
    else
    {
      op.cps <- c(op.cps, max(which((criterion) == TRUE)))
    }
    if (op.cps == Q)
    {
      warning("The number of changepoints identified is Q, it is advised to increase Q to make sure changepoints have not been missed.")
    }
    if (op.cps > 0)
    {
      cpts = sort(cpt[1, 1:op.cps])
    }
  }
  return(cpts)
}

#' Stationarity checker and divider
#'
#' Finds all stationar periods of data. Maximal number of intervals is nchanges
#' @author Harutyun Khachatryan
#' @param itseries time series like vector
#' @param iobs length of window around ind index
#' @param xhint logical TRUE/FALSE indicating if data is highly unstable (default value is FALSE)
#' @param xper logical TRUE/FALSE indicating if data is periodic (default value is FALSE)
#' @param nchanges number of possible changes in data (default value is NA,
#' in that case nchanges is equal to one third of number of days in data)
#' @export
StatLen2 <- function(itseries, iobs, xhint = FALSE, xper = FALSE, nchanges = NA)
{
  intseries <- itseries
  l3dstat <- 0
  trend_det <- 0
  nobsi <- length(itseries)
  stati2 <- list(strong = c(1, nobsi), weak = c(1, nobsi))
  obsi <- iobs
  if ((nobsi - 2 * obsi) >= 0)
  {
    rnch <- round(nobsi / (3 * obsi))
  } else
  {
    rnch <- 0
  }
  if ( (!is.na(nchanges)) & (obsi >= 24) )
  {
    rnch <-  nchanges
  }
  if (rnch > 0)
  {
    error.message <- try(lchanges <- MultipleChange2(as.vector(intseries), Q = rnch), silent = TRUE)
    if ((length(lchanges) < rnch) | (!xper))
    {
      error.message <- try(xchange <- SingleChange2(as.vector(intseries)), silent = TRUE)
    } else
    {
      xchange <- NULL
    }
    if (!xper)
    {
      if (nobsi > 3 * obsi)
      {
        trend_det <- STest2(intseries[(nobsi - 2 * obsi):nobsi]) * STest2(intseries[(nobsi - obsi):nobsi])
        l3dstat <- MVTest2(intseries, (nobsi - 2 * obsi), 2 * obsi) + ATest2(intseries, (nobsi - 2 * obsi), 2 * obsi) - trend_det
      }
      if ((nobsi > 3 * obsi) & (l3dstat > 0))
      {
        error.message <- try(lchanges0 <- (MultipleChange2(as.vector(intseries[(nobsi - 2 * obsi - 1):nobsi]), Q = 3, robust = TRUE) + (nobsi - 2 * obsi - 2)), silent = TRUE)
      }
      if ( (exists("lchanges")) & (exists("lchanges0")) )
      {
        if ( (length(lchanges) > 0) & (length(lchanges0) > 0) )
        {
          dchanges <- which(sapply(sapply(1:length(lchanges0), function(w) which(abs(lchanges0[w] - lchanges) < 3)), length) > 0)
          if (length(dchanges) > 0) {lchanges0 <- lchanges0[-dchanges]}
        }
        if (length(lchanges0) > 0) {lchanges <- sort(c(lchanges, lchanges0))}
      }
      if ( (!exists("lchanges")) & (exists("lchanges0")) )
      {
        lchanges <- lchanges0
      }
      if (exists("lchanges"))
      {
        if ( (length(lchanges) > 0) & (length(xchange) > 0) & (min(abs(lchanges - xchange), na.rm = TRUE) >= 3) )
        {
          lchanges <- sort(c(lchanges, xchange))
        } else
        {
          if ( (length(lchanges) == 0) & (length(xchange) > 0) )
          {
            lchanges <- xchange
          }
        }
      }
    } else
    {
      if ( (length(lchanges) > 0) & (length(xchange) > 0) )
      {
        if (!any(lchanges == xchange))
        {
          lchanges <- sort(c(lchanges, xchange))
        }
      }
      if ( (length(lchanges) == 0) & (length(xchange) > 0) )
      {
        lchanges <- xchange
      }
    }
  }
  if (exists("lchanges"))
  {
    if (length(lchanges) > 0)
    {
      if ((trend_det > 0) & (any(lchanges > (nobsi - obsi))))
      {
        lchanges <- lchanges[-which(lchanges > (nobsi - obsi))]
      }
      ltest <- sapply(lchanges, function(w) {ATest2(intseries, w, obsi) + MVTest2(intseries, w, obsi) + KTest2(intseries, w, obsi, 2)})
      if (!xhint)
      {
        kltest <- sapply(lchanges, function(w) KTest2(intseries, w, obsi, 1))
      }
      if (xhint)
      {
        kltest <- sapply(lchanges, function(w) KTest2(intseries, w, obsi, 0))
      }
      if ( (!xhint) & (any(ltest > 0.19)) )
      {
        stati2$strong <- c(1, lchanges[which(ltest > 0.19)], nobsi)
      }
      if ( (xhint) & (any(ltest > 0.06)) )
      {
        stati2$strong <- c(1, lchanges[which(ltest > 0.06)], nobsi)
      }
      if (any(kltest > 0.06))
      {
        stati2$weak <- c(1, lchanges[which(kltest > 0.06)], nobsi)
      }
    }
  }
  stati2$weak <- as.integer(stati2$weak)
  stati2$strong <- as.integer(stati2$strong)
  return(stati2)
}

#' Bad data finder
#'
#' Finds locations of bad days in data according to stationarity change points
#' @author Harutyun Khachatryan
#' @param intseries time series like vector
#' @param norm_dates timeseries observation times given as vector of POSIXct
#' @param real_dates real observation times given as vector of POSIXct
#' @param stat_changes indeces of stationarity changes in time series
#' @param xhint logical TRUE/FALSE indicating if data is highly unstable (default value is FALSE)
#' @export
BadInd2 <- function(intseries, norm_dates, real_dates, stat_changes, xhint = FALSE)
{
  nobsi <- length(intseries)
  delta_t <- as.numeric(norm_dates[2]) - as.numeric(norm_dates[1])
  obsi <- ceiling(24 * 60 * 60 / delta_t)
  schanges <- stat_changes
  erase_inf <- list(indices = NULL, stat = schanges)
  if (length(schanges) > 0)
  {
    binds <- which(diff(schanges) < obsi)
    if (length(binds) > 0)
    {
      erase_inds <- mat.or.vec(0,0)
      rerase_inds <- mat.or.vec(0,0)
      chmark <- mat.or.vec(0,length(binds))
      for (i in 1:length(binds))
      {
        if ((schanges[binds[i]] + 1 + obsi) < (nobsi - obsi))
        {
          chts <- intseries[-seq(schanges[binds[i]] + 1, schanges[binds[i]] + 1 + obsi)]
        }
        else
        {
          chts <- intseries[-seq(schanges[binds[i]] + 1, schanges[binds[i] + 1] - 1)]
        }
        if (KTest2(chts, schanges[binds[i]], obsi, 2) == 0) {chmark[i] <- ATest2(chts, schanges[binds[i]], obsi) + MVTest2(chts, schanges[binds[i]], obsi)} else {chmark[i] <- 0.2}
        if ((!xhint) & (chmark[i] == 0.1)) {chmark[i] <- 0}
        if (chmark[i] < 0.06)
        {
          if ((schanges[binds[i]] + obsi) <= nobsi)
          {
            erase_inds <- append(erase_inds, seq(schanges[binds[i]], (schanges[binds[i]] + obsi)))
          } else
          {
            erase_inds <- append(erase_inds, seq(schanges[binds[i]], nobsi))
          }
          remdate <- real_dates[max(which(real_dates <= norm_dates[schanges[binds[i]]]), na.rm = TRUE)]
          rerase_inds <-  append(rerase_inds, which(((real_dates <= (remdate + 24 * 60 * 60)) & (real_dates > remdate))))
        }
      }
      erase_inds <- unique(erase_inds)
      rerase_inds <- unique(rerase_inds)
      if (any(chmark < 0.06)) {schanges <- schanges[-c(binds[which(chmark < 0.06)], (binds[which(chmark < 0.06)] + 1))]}
      erase_stat <- intersect(erase_inds, schanges)
      if (length(erase_stat) > 0)
      {
        schanges <- setdiff(schanges, erase_stat)# damn assymetric setdiff!!
      }
      erase_inf$indices <- as.integer(rerase_inds)
      erase_inf$stat <- schanges
    }
  }
  return(erase_inf)
}
