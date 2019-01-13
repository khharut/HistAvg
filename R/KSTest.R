
#' Estimates Kolmogorov-Smirnov's test statistic lambda, wrapper
#'
#' Extracts from stats ks.test method main statistic lambda, works on any input data
#' @author Harutyun Khachatryan
#' @param itseries time series like vector of values
#' @param ind index around which test is performed
#' @param lobs length of window around ind index
#' @return Kolmogorov-Smirnov's test lambda
#' @export
#' @importFrom stats ks.test
KSWrapper2 <- function(itseries, ind, lobs)
{
  lambda <- 0
  ktstart <- ind - lobs + 1
  ktend <- ind + lobs
  nobsi <- length(itseries)
  if ( (ktend <= nobsi) & (ktstart >= 1) )
  {
    lcoef <- sqrt(lobs / 2)
    try(test_res <- suppressWarnings(ks.test(itseries[ktstart:ind], itseries[(ind + 1):ktend])), silent = TRUE)
    if (exists("test_res"))
    {
      lambda <- lcoef * test_res$statistic
      lambda <- unname(lambda)
    }
  }
  return(lambda)
}

#' Estimates Kolmogorov-Smirnov's test statistic lambda
#'
#' Gives a two-dimensional matrix of indeces and lambdas
#' @author Harutyun Khachatryan
#' @param itseries time series like vector of values
#' @param lobs length of window around ind index
#' @return two-dimensional matrix of indeces and lambdas
#' @export
KSValues2 <- function(itseries, lobs)
{
  intseries <- itseries
  nobsi <- length(itseries)
  obsi <- lobs
  ksvalues <- cbind(c(0, nobsi), c(0, 0))
  if ((nobsi - 2 * obsi) >= 0)
  {
    kindeces <- seq((obsi + 1), nobsi, obsi)
    if (kindeces[length(kindeces)] == nobsi)
    {
      kindeces <- kindeces[-length(kindeces)]
    }
    else
    {
      if (kindeces[length(kindeces)] > (nobsi - obsi))
      {
        kindeces[length(kindeces)] <- (nobsi - obsi)
      } else
      {
        kindeces <- c(kindeces, (nobsi - obsi))
      }
    }
    if (obsi >= 24)
    {
      lambdas <- sapply(kindeces, function(w) KSWrapper2(intseries, w, obsi))
    } else {lambdas <- rep(0, length(kindeces))}
    ksvalues <- cbind(kindeces, lambdas, deparse.level = 0)
  }
  if (!is.matrix(ksvalues))
  {
    ksvalues <- as.matrix(ksvalues)
    ksvalues <- t(ksvalues)
  }
  return(ksvalues)
}

#' Performs Kolmogorov-Smirnov's test
#'
#' Tests if ind index is a point of stationarity change by mean of Kolmogorov-Smirnov test
#' @author Harutyun Khachatryan
#' @param itseries time series like vector of values
#' @param ind index around which test is performed
#' @param lobs length of window around ind index
#' @param robust integer value with level 0, 1, 2 (default 1)
#' a) when robust is 0 then estimated KS lambda should be greater or equal than exp(1) for detecting ind as change point
#' b) when robust is 1 then estimated KS lambda should be greater or equal than 4.2 for detecting ind as change point
#' c) when robust is 2 then estimated KS lambda should be greater or equal than 7.0 for detecting ind as change point
#' @return 0 if no stationarity change is detected 0.1 otherwise
#' @export
KTest2 <- function(itseries, ind, lobs, robust = 1)
{
  cpvalue <- 0
  nobsi <- length(itseries)
  ktest <- 0
  if ((ind > lobs) & (ind < (nobsi - lobs + 2)))
  {
    ktest <- KSWrapper2(itseries, ind, lobs)
  } else
  {
    if (nobsi >= (2 * lobs + 2))
    {
      if (ind <= lobs)
      {
        ktest <- KSWrapper2(itseries, lobs + 1, lobs)
      }
      if (ind >= (nobsi - lobs + 2))
      {
        ktest <- KSWrapper2(itseries, nobsi - lobs + 1, lobs)
      }
    }
  }
  if ((robust == 0) & (ktest >= exp(1))) {cpvalue <- 0.1}
  if ((robust == 1) & (ktest >= 4.2)) {cpvalue <- 0.1}
  if ((robust == 2) & (ktest >= 7)) {cpvalue <- 0.15}## exp(2) should be changed to 7
  return(cpvalue)
}

#' My own redifintion of round function
#'
#' Redifinition of round function in cases when even number +0.5 rounded downstair
#' @author Harutyun Khachatryan
#' @param x numeric value to be rounded
#' @return rounded value xr of original value x
#' @export
HRound <- function(x)
{
  xr1 <- round(x)
  xr2 <- round(x + 1)
  if ((is.finite(xr1)) & (is.finite(xr2)))
  {
    if ((xr2 - xr1) > 1.5)
    {
      xr <- xr2 - 1
    } else
    {
      xr <- xr1
    }
    return(xr)
  } else
  {
    return(x)
  }
}

#' Extra information from Kolmogorov-Smirnov lambda
#'
#' Finds different information from KS lambdas like: number of change points, weekly pattern, stationarity
#'
#' @author Harutyun Khachatryan
#' @param  lambdas 2-dimensional matrix of KS lambdas and indeces
#' @return a list of different parameters:
#' 1) a field named "weekchange" conatins indeces of stationarity change points affected by week pattern existance in data
#' 2) a field named "type" has two possible values "low" and "high" indicating type of week pattern affected changes
#' 3) a field named "nchange" number of possible stationarity changes in data
#' 4) a field named "weekly7" logical indicating if just 7 day pattern exists
#' 5) a field named "weekly52" logical indicating if just 5-2 (working/non-working) days pattern exists
#' 6) a field named "highly_unstable" logical showing if data is higly unstable, can have increasing or decreasing trend
#' @export
KSInfo2 <- function(lambdas)
{
  weekper <- 0
  ee <- exp(1)
  if (!is.matrix(lambdas))
  {
    lambdas <- as.matrix(lambdas)
    lambdas <- t(lambdas)
  }
  ksinfo <- list(weekchange = NULL, type = 0, nchange = 0, weekly7 = c(FALSE, FALSE), weekly52 = c(FALSE, FALSE), highly_unstable = FALSE)
  if (any(lambdas[, 2] > 4)) {ee <- 2.5}
  change_num <- sum(lambdas[, 2] >= ee)
  ksinfo$nchange <- change_num
  ks_len <- length(lambdas[, 2])
  if (length(lambdas[ ,1]) > 1)
  {
    obsi <- lambdas[2, 1] - lambdas[1, 1]
    nobsi <- max(lambdas[, 1], na.rm = TRUE) + obsi
    weeksn <- nobsi / (7 * obsi)
    weeksnt <- HRound(weeksn)
  } else
  {
    obsi <- 0
    nobsi <- 0
    weeksn <- 0
    weeksnt <- 0
  }
  if ((change_num == ks_len) & (change_num > 0)) {ksinfo$highly_unstable <- TRUE}
  if ( (weeksn > 2) & (change_num > 0) )
  {
    if ( (lambdas[ks_len, 1] - lambdas[(ks_len - 1), 1]) < obsi )
    {
      high_ks <- lambdas
      high_ks <- high_ks[-ks_len, ]
    } else {high_ks <- lambdas}
    weekper <- NRel2(high_ks[, 2], 7)[2]
    if (length(high_ks[, 1]) > (7 * weeksnt))
    {
      high_ks <- high_ks[(length(high_ks[, 1]) - 7 * weeksnt + 1):length(high_ks[, 1]),]
    }
    low_ks <- high_ks
    if ( length(which( high_ks[, 2] >= ee )) > 1 )
    {
      high_ks <- high_ks[which( high_ks[, 2] >= ee ),]
      ksorder <- order(high_ks[, 2], decreasing = TRUE)
      high_ks <- high_ks[ksorder, ]
      if ( length(high_ks[, 2]) >= weeksnt )
      {
        pattern7d <- sort(high_ks[1:weeksnt, 1])
        kspattern7d <- diff(pattern7d)
        kspattern7d <- kspattern7d / obsi
        ksinfo$weekly7[1] <- all(kspattern7d == 7)
      }
      if ( length(high_ks[, 2]) >= (2 * weeksnt) )
      {
        pattern52d <- sort(high_ks[1:(2 * weeksnt), 1])
        kspattern52d <- diff(pattern52d)
        kspattern52d <- kspattern52d / obsi
        kspattern52d <- unique(kspattern52d)
        if (length(kspattern52d) == 2)
        {
          ksinfo$weekly52[1] <- ((kspattern52d[1] + kspattern52d[2]) == 7)
        }
      }
    }
    if ( length(which( low_ks[, 2] < ee )) > 1 )
    {
      low_ks <- low_ks[which( low_ks[, 2] < ee ),]
      lksorder <- order(low_ks[, 2])
      low_ks <- low_ks[lksorder, ]
      if ( length(low_ks[, 2]) >= weeksnt )
      {
        lpattern7d <- sort(low_ks[1:weeksnt, 1])
        lkspattern7d <- diff(lpattern7d)
        lkspattern7d <- lkspattern7d / obsi
        ksinfo$weekly7[2] <- all(lkspattern7d == 7)
      }
      if ( length(low_ks[, 2]) >= (2 * weeksnt) )
      {
        lpattern52d <- sort(low_ks[1:(2 * weeksnt), 1])
        lkspattern52d <- diff(lpattern52d)
        lkspattern52d <- lkspattern52d / obsi
        lkspattern52d <- unique(lkspattern52d)
        if (length(lkspattern52d) == 2)
        {
          ksinfo$weekly52[2] <- ((lkspattern52d[1] + lkspattern52d[2]) == 7)
        }
      }
    }
  }
  if ((ksinfo$weekly7[1]) | (ksinfo$weekly52[1]))
  {
    ksinfo$type <- 3
  }
  else
  {
    if ((ksinfo$weekly7[2]) | (ksinfo$weekly52[2]))
    {
      ksinfo$type <- 2
    }
  }
  if (ksinfo$weekly52[1]) {ksinfo$weekchange <- pattern52d}
  if ((ksinfo$weekly7[1]) & (length(ksinfo$weekchange) == 0)) {ksinfo$weekchange <- pattern7d}
  ksinfo$weekly52 = ksinfo$weekly52[1] | ksinfo$weekly52[2]
  ksinfo$weekly7 = ksinfo$weekly7[1] | ksinfo$weekly7[2]
  if ((ksinfo$type == 0) & (weekper > 0.75) & (change_num >= weeksn))
  {
    ksinfo$type <- 1
  }
  if ((ksinfo$type > 0) & (weeksn >= 3))
  {
    if ( (lambdas[ks_len, 1] - lambdas[(ks_len - 1), 1]) < obsi )
    {
      klambdas <- lambdas
      klambdas <- klambdas[-ks_len, ]
    } else {klambdas <- lambdas}
    ktseries <- ts(klambdas[, 2], frequency = 7)
    kses <- stl(ktseries, s.window = "per")$time.series[1:7, "seasonal"]
    kses_order <- order(kses, decreasing = TRUE)
    apattern7d <- kses_order[1] %% 7
    apattern52d <- kses_order[2] %% 7
    if ((abs(apattern7d - apattern52d) <= 2) & (ksinfo$weekly7))
    {
      ksinfo$weekly7 <- FALSE
      ksinfo$weekly52 <- TRUE
      length(ksinfo$weekchange) <- 0
    }
    apattern7d <- kses_order[1]
    apattern52d <- kses_order[2]
    apattern7d <- seq(apattern7d, (7 * weeksn), 7)
    apattern52d <- seq(apattern52d, (7 * weeksn), 7)
    apattern52d <- sort(c(apattern52d, apattern7d))
    klindeces <- lambdas[, 1]
    apattern7d <- klindeces[apattern7d]
    apattern52d <- klindeces[apattern52d]
    if ( (ksinfo$weekly52) & (length(ksinfo$weekchange) == 0) ) {ksinfo$weekchange <- apattern52d}
    if ( (ksinfo$weekly7) & (length(ksinfo$weekchange) == 0) ) {ksinfo$weekchange <- apattern7d}
    if (length(ksinfo$weekchange) == 0)
    {
      ksinfo$weekchange <- apattern52d
    }
  }
  if ((ksinfo$type != 0) & (ksinfo$nchange %% 2))
  {
    ksinfo$nchange <- ksinfo$nchange + 1
  }
  if ((length(ksinfo$weekchange) >= ksinfo$nchange) & (length(ksinfo$weekchange) > 0) & (ksinfo$type == 3))
  {
    ksinfo$nchange <- 0
  }
  return(ksinfo)
}
