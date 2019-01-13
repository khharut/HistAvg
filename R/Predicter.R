
#' Creates linear prediction model for given data
#'
#' Constructs linear model containing season and/or trend
#' @author Harutyun Khachatryan
#' @param data timseries like data to be modeled
#' @param lobs number of observations per day
#' @param steps number of data points to be predicted
#' @param trend_only logical TRUE/FALSE indicates if seasonal part should be included
#' @param sensitivity shows confidence interval has three values "low", "medium" and "high" (default is "low")
#' @return prediction based on linear model and confidence intervals along with model reliability
#' @export
LinearModel2 <- function(data, lobs, steps = lobs, trend_only = FALSE, sensitivity = "low")
{
  R_square <- 0
  if ((sensitivity != "low") & (sensitivity != "medium") & (sensitivity != "high") )
  {
    sensitivity <- "low"
  }
  if (sensitivity == "low")
  {
    pvalue <- 99
  }
  if (sensitivity == "medium")
  {
    pvalue <- 95
  }
  if (sensitivity == "high")
  {
    pvalue <- 90
  }
  levels <- c(50, pvalue)
  nobsi <- length(data)
  trend <- seq(1, nobsi + steps)
  if (!trend_only)
  {
    season <- rep(seq(1, lobs), length.out = nobsi + steps)
    season <- as.factor(season)
    newdata <- data.frame(trend, season)
    newdata <- newdata[(nobsi + 1):(nobsi + steps), ]
    length(season) <- nobsi
    length(trend) <- nobsi
    linear_model <- lm(data ~ season + trend)
  } else
  {
    newdata <- list()
    newdata$trend <- trend[(nobsi + 1):(nobsi + steps)]
    length(trend) <- nobsi
    linear_model <- lm(data ~ trend)
  }
  uncov_var <- sum(linear_model$residuals ^ 2)
  cov_var <- sum(((linear_model$fitted.values) - mean(linear_model$fitted.values)) ^ 2)
  if ((cov_var + uncov_var) > 0)
  {
    R_square <- cov_var / (cov_var + uncov_var)
    R_square <- abs(1 - ( (1 - R_square) * (nobsi - 1) / (nobsi - lobs - 1) ))
  } else
  {
    R_square <- 1.0
  }
  linear_prediction50 <- predict(linear_model, newdata = newdata, level = levels[1] / 100, interval = "prediction", se.fit = TRUE)$fit
  linear_prediction99 <- predict(linear_model, newdata = newdata, level = levels[2] / 100, interval = "prediction", se.fit = TRUE)$fit
  data_prediction <- list(mean = linear_prediction50[, 1], lower = cbind(linear_prediction50[, 2], linear_prediction99[, 2]),
                         upper = cbind(linear_prediction50[, 3], linear_prediction99[, 3]), reliability = sqrt(R_square))
  return(data_prediction)
}

#' Estimates prediction based dynamic thresholds
#'
#' Develops dynamics thresholds using simple lienear model with trend and season
#' @author Harutyun Khachatryan
#' @param norm_data timeseries like observed with fixed time step values
#' @param norm_times timeseries observation times given as vector of POSIXct
#' @param is_positive logical TRUE/FALSE indicating that data non-negative or not
#' @param trend_only logical TRUE/FALSE indictaing if only trend will be predicted
#' @param psense shows confidence interval has three values "low", "medium" and "high" (default is "low")
#' @return prediction dates along with lower, upper limits for data and reliability
#' @export
#' @importFrom forecast BoxCox.lambda
Predicter2 <- function(norm_data, norm_times, is_positive = TRUE, trend_only = FALSE, psense = "low")
{
  last_day <- trunc(norm_times[length(norm_times)], "day")
  delta_t <- as.numeric(norm_times[2]) - as.numeric(norm_times[1])
  delta_t <- delta_t * 1000
  obsi = ceiling( (24 * 60 * 60 * 1000) / delta_t)
  pred_dates <- seq(from = norm_times[length(norm_times)] + (delta_t / 1000),
                   to = last_day + 2 * 24 * 60 * 60, by = paste(as.character(delta_t / 1000), "secs", sep = " ", collapse = NULL)) - 24 * 60 * 60
  linear_prediction <- list(thresholds = cbind(pred_dates, rep(0, length(pred_dates)), rep(0, length(pred_dates))), reliability = 0)
  if (pred_dates[length(pred_dates)] == (last_day + 24 * 60 * 60) ) {pred_dates <- pred_dates[-length(pred_dates)]}
  mlambda <- BoxCox.lambda(norm_data)
  x_forecast <- LinearModel2(YeoJohn2(norm_data, mlambda), obsi, steps = length(pred_dates), trend_only = trend_only, sensitivity = psense)
  x_forecast$mean <- YeoJohnInv2(x_forecast$mean, mlambda)
  x_forecast$upper[,2] <- YeoJohnInv2(x_forecast$upper[, 2], mlambda)
  x_forecast$lower[,2] <- YeoJohnInv2(x_forecast$lower[, 2], mlambda)
  x_forecast$lower[,1] <- YeoJohnInv2(x_forecast$lower[, 1], mlambda)
  x_forecast$upper[,1] <- YeoJohnInv2(x_forecast$upper[, 1], mlambda)
  if ( (any(is.na(x_forecast$upper[, 2]))) | (any(is.na(x_forecast$lower[, 2]))) )
  {
    if ( (length(which(!is.na(x_forecast$upper[, 2]))) > 1) & (length(which(!is.na(x_forecast$lower[, 2]))) > 1) )
    {
      divx <- as.vector(x_forecast$upper[, 2] - x_forecast$lower[, 2])
      if (all(!is.na(x_forecast$mean)))
      {
        divm <- x_forecast$mean
      } else
      {
        if (all(!is.na(x_forecast$lower[, 1]))) {divm <- x_forecast$lower[, 1]}
      }
      if (!exists("divm"))
      {
        divm <- approx(x = as.numeric(names(x_forecast$mean)),y = x_forecast$mean, xout = as.numeric(names(x_forecast$mean)), method = "linear", rule = 2)$y
      }
      div_lm <- lm(divx ~ divm)
      divx_est <- div_lm$coefficients[1] + divm * div_lm$coefficients[2]
      if (any(is.na(divx)))
      {
        divx[is.na(divx)] <- divx_est[is.na(divx)]
      }
      if (any(is.na(x_forecast$upper[, 2])))
      {
        x_forecast$upper[is.na(x_forecast$upper[, 2])] <- x_forecast$lower[is.na(x_forecast$upper[, 2]), 2] + divx[is.na(x_forecast$upper[, 2])]
      }
      if (any(is.na(x_forecast$lower[, 2])))
      {
        x_forecast$lower[is.na(x_forecast$lower[, 2])] <- x_forecast$upper[is.na(x_forecast$lower[, 2]), 2] - divx[is.na(x_forecast$lower[, 2])]
      }
    } else
    {
      if ((length(which(!is.na(x_forecast$lower[, 2]))) > 1))
      {
        div_lm_low <- lm(x_forecast$lower[, 2] ~ x_forecast$lower[, 1])
        upper_est <- (x_forecast$upper[, 1] - div_lm_low$coefficients[1]) / div_lm_low$coefficients[2]
        if (any(is.na(x_forecast$upper[, 2])))
        {
          x_forecast$upper[is.na(x_forecast$upper[, 2]), 2] <- upper_est[is.na(x_forecast$upper[, 2])]
        }
      }
      if ((length(which(!is.na(x_forecast$upper[, 2]))) > 1))
      {
        div_lm_up <- lm(x_forecast$upper[, 2] ~ x_forecast$upper[, 1])
        lower_est <- (x_forecast$lower[, 1] - div_lm_up$coefficients[1]) / div_lm_up$coefficients[2]
        if (any(is.na(x_forecast$lower[, 2])))
        {
          x_forecast$lower[is.na(x_forecast$lower[, 2]), 2] <- lower_est[is.na(x_forecast$lower[, 2])]
        }
      }
    }
  }
  x_predict_mean <- as.vector(x_forecast$mean)
  x_predict_ulevel99 <- as.vector(x_forecast$upper[, 2])
  x_predict_llevel99 <- as.vector(x_forecast$lower[, 2])
  x_forecast$upper <- x_forecast$upper[, -1]
  x_forecast$lower <- x_forecast$lower[, -1]
  if (is_positive)
  {
    if (length(which(x_predict_mean < 0)) > 0)
    {
      x_predict_mean[which(x_predict_mean < 0)] <- 0
    }
    if (length(which(x_predict_ulevel99 < 0)) > 0)
    {
      x_predict_ulevel99[which(x_predict_ulevel99 < 0)] <- 0
    }
    if (length(which(x_predict_llevel99 < 0)) > 0)
    {
      x_predict_llevel99[which(x_predict_llevel99 < 0)] <- 0
    }
  }
  linear_prediction <- list(thresholds = cbind(pred_dates, x_predict_ulevel99, x_predict_llevel99), reliability = x_forecast$reliability)
  if (trend_only)
  {
    linear_prediction$thresholds <- linear_prediction$thresholds[seq((length(pred_dates) - obsi + 1), length(pred_dates)), ]
  }
  return(linear_prediction)
}

#' Performs Yeo-Johnson power transformation on a vector
#'
#' Does Yeo-Johnson transformation on a vector with parameter lambda
#' @author Harutyun Khachatryan
#' @param x vector or timeseries to be transformed
#' @param lambda transformation parameter obtained from BoxCox.lambda method of forecast package
#' @return transformed values of input vector
#' @export
YeoJohn2 <- function(x, lambda = 0)
{
  out <- x
  if ( (any(x >= 0)) & (lambda != 0) ) {out[x >= 0] <- ((((out[x >= 0] + 1) ^ (lambda)) - 1) / lambda)}
  if ( (any(x >= 0)) & (lambda == 0) ) {out[x >= 0] <- log(out[x >= 0] + 1)}
  if ( (any(x < 0)) & (lambda != 2) ) {out[x < 0] <- -((((1 - out[x < 0]) ^ (2 - lambda)) - 1) / (2 - lambda))}
  if ( (any(x < 0)) & (lambda == 2) ) {out[x < 0] <- -log(1 - out[x < 0])}
  return(out)
}

#' Performs Yeo-Johnson inverse power transformation on a vector
#'
#' Does invserse (i.e. back) Yeo-Johnson transformation on a vector with parameter lambda
#' @author Harutyun Khachatryan
#' @param x vector or timeseries to be transformed
#' @param lambda transformation parameter obtained from BoxCox.lambda method of forecast package
#' @return inverse transformed values of input vector
#' @export
YeoJohnInv2 <- function(x, lambda = 0)
{
  out <- x
  if ( (any(x >= 0)) & (lambda != 0) ) {out[x >= 0] <- ((out[x >= 0] * lambda + 1) ^ (1 / lambda)) - 1}
  if ( (any(x >= 0)) & (lambda == 0) ) {out[x >= 0] <- exp(out[x >= 0]) - 1}
  if ( (any(x < 0)) & (lambda != 2) ) {out[x < 0] <- 1 - ((1 - out[x < 0] * (2 - lambda)) ^ (1 / (2 - lambda)))}
  if ( (any(x < 0)) & (lambda == 2) ) {out[x < 0] <- 1 - exp(-out[x < 0])}
  return(out)
}

#' Predicts mean value
#'
#' Performs mean value prediction and return confidence interval for it
#' @author Harutyun Khachatryan
#' @param x values to be predict
#' @param sensitivity shows confidence interval has three values "low", "medium" and "high" (default is "low")
StepY2 <- function(x, sensitivity = "low")
{
  if ((sensitivity != "low") & (sensitivity != "medium") & (sensitivity != "high") )
  {
    sensitivity <- "low"
  }
  if (sensitivity == "low")
  {
    pvalue <- 99
  }
  if (sensitivity == "medium")
  {
    pvalue <- 95
  }
  if (sensitivity == "high")
  {
    pvalue <- 90
  }
  n <- length(x)
  meanx <- mean(x, na.rm = TRUE)
  f <- meanx
  lower <- upper <- NA
  s <- sd(x, na.rm = TRUE)
  tfrac <- qt(0.5 - pvalue / 200, n - 1)
  w <- -tfrac * s * sqrt( 1 + 1 / n)
  lower <- f - w
  upper <- f + w
  return(c(lower, upper))
}
