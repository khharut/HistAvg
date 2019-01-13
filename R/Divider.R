
#' Projects all data to one day and divides data into groups
#'
#' Takes timestamps and values and creates groups according to projected timestamp position in day
#' @author Harutyun Khachatryan
#' @param time_stamps one dimensional array of timestamp
#' @param values one dimensional array of values measured at times 'time_stamps', i.e. they should have the same length
#' @param delta_min time aggregation parameter given in minutes
#' @param ch_ind integer value of index at which concept changes last occur
#' @param overlap logical TRUE/FALSE, (default TRUE) indicating if derived groups will be overlapping or not
#' @param loop logical TRUE/FALSE, (default TRUE) indicating if last and first groups will be mirrored or not
#' @return list of vectors representing values in groups, name of each list is timestamps around which values appeared in one day projection
#' @note Monitis
#' @keywords time series
#' @export
#' @importFrom gmp as.bigz
DivGroups2 <- function(time_stamps, values, delta_min, ch_ind = NA, overlap = TRUE, loop = TRUE)
{
  day_mili <- 24 * 60 * 60 * 1000
  delta_mili <- delta_min * 60000
  ch_ind <- as.numeric(ch_ind)
  if (!is.na(ch_ind))
  {
    if (ch_ind > 0)
    {
      time_stamps <- time_stamps[(ch_ind + 1):length(time_stamps)]
      values <- values[(ch_ind + 1):length(values)]
    }
  }
  x_groups <- list()
  last_day <- day_mili * trunc(max(time_stamps, na.rm = TRUE) / day_mili)
  if (delta_mili <= (day_mili / 4))
  {
    oneday <- time_stamps %% day_mili
    if (overlap)
    {
      evens <- trunc(oneday / delta_mili) + 1
      evens <- 2 * evens
      odds <- trunc((oneday + (0.5 * delta_mili)) / delta_mili)
      odds <- 2 * odds + 1
      all_ind <- c(evens, odds)
      ngroup <- ceiling(2 * day_mili / delta_mili) + 1
      time_shift <- (ngroup - 1) * (0.5 * delta_mili) - day_mili
      nod_ind <- seq(1, ngroup)
      split_nodes <- 0.5 * delta_mili * (nod_ind - 1)
      x_groups[as.character(nod_ind)] <- NA
      x_groups[sort(unique(all_ind))] <- split(c(values, values), all_ind)
      if (time_shift > 0)
      {
        end_group <- which(oneday >= (day_mili - 0.5 * delta_mili))
        x_groups[[ngroup]] <- values[end_group]
        split_nodes[ngroup] <- day_mili
      }
    } else
    {
      all_ind <- trunc((oneday + (0.5 * delta_mili)) / delta_mili) + 1
      ngroup <- ceiling(day_mili / delta_mili) + 1
      time_shift <- (ngroup - 1) * (delta_mili) - day_mili
      nod_ind <- seq(1, ngroup)
      split_nodes <- delta_mili * (nod_ind - 1)
      x_groups[as.character(nod_ind)] <- NA
      x_groups[sort(unique(all_ind))] <- split(values, all_ind)
      if (time_shift > 0)
      {
        end_group <- which(oneday >= (day_mili - delta_mili))
        x_groups[[ngroup]] <- values[end_group]
        split_nodes[ngroup] <- day_mili
      }
    }
    if (loop)
    {
      if ((!all(is.na(x_groups[[1]]))) & (!all(is.na(x_groups[[length(x_groups)]]))))
      {
        x_groups[[1]] <- c(x_groups[[1]], x_groups[[length(x_groups)]])
        x_groups[[length(x_groups)]] <- x_groups[[1]]
      } else
      {
        if (!all(is.na(x_groups[[1]]))) {x_groups[[length(x_groups)]] <- x_groups[[1]]}
        if (!all(is.na(x_groups[[length(x_groups)]]))) {x_groups[[1]] <- x_groups[[length(x_groups)]]}
        if ((all(is.na(x_groups[[1]]))) & (all(is.na(x_groups[[length(x_groups)]]))))
        {
        	length(x_groups[[1]]) <- 1
        	length(x_groups[[length(x_groups)]]) <- 1
        }
      }
    }
    split_nodes <- split_nodes + last_day
    names(x_groups) <- as.character(as.bigz(split_nodes))
  } else
  {
    x_groups[as.character(gmp::as.bigz(last_day))] <- NA
  }
  return(x_groups)
}

#' Projects timestamps to one day, i.e. last appeared day in given vector of timestamps
#'
#' Takes original timestamps and projects them all to last day appeared in timestamps
#' @author Harutyun Khachatryan
#' @param time_stamps original timestamps
#' @param ch_ind integer value of index at which concept changes last occur
#' @return one day projected timestamps, i.e. some values of timestamps is changed
#' @export
OneDayProj2 <- function(time_stamps, ch_ind = NA)
{
  day_mili <- 24 * 60 * 60 *1000
  if (!is.na(ch_ind))
  {
    if (ch_ind > 0)
    {
      time_stamps <- time_stamps[(ch_ind + 1):length(time_stamps)]
    }
  }
  oneday <- time_stamps %% day_mili
  return(oneday)
}

#' Calculates quantiles for historical data
#'
#' Gives quantiles of historical data for predefined quantiles hprobs
#' @author Harutyun Khachatryan
#' @param x historical data values
#' @param hprobs probs of quantiles needed for calculation, for example c(0.01,0.99)
#' @param limits limits on values within which most of historical values should lay but for some caes it can goes out
#' @param nout integer value shows how many historical values should go out limits that it can be treated as tendency
#' @return lower and upper bounds for historical data
#' @export
HistQuantile2 <- function(x, hprobs, limits, nout)
{
  qdown <- 0
  qup <- 0
  x.sd <- sd(x, na.rm = TRUE)
  if ( (x.sd > 0) & (!is.na(x.sd)) )
  {
    hquant <- as.vector(quantile(x, probs = hprobs, na.rm = TRUE))
    qupper <- max(hquant, na.rm = TRUE)
    qlower <- min(hquant, na.rm = TRUE)
    quplow <- c(max(x[which(x <= qupper)]), min(x[which(x >= qlower)]))
    qup <- max(quplow, na.rm = TRUE)
    qdown <- min(quplow, na.rm = TRUE)
    if ( (qup > limits[2]) & (length(which(((x > limits[2]) & (x < qup)))) < nout) )
    {
      qup <- max(x[x <= limits[2]])
    }
    if ((qdown < limits[1]) & (length(which(((x < limits[1]) & (x > qdown)))) < nout) )
    {
      qdown <- min(x[x >= limits[1]])
    }
  } else
  {
    qup <- mean(x, na.rm = TRUE)
    qdown <- mean(x, na.rm = TRUE)
  }
  return(c(qdown, qup))
}

#' Gives an array of outliers
#'
#' Finds outliers from given limits for historical data and represents every outlier by it's distance from limits
#' @author Harutyun Khachatryan
#' @param div_groups one day projected groups of data in list of vectors
#' @param limits a list of vectors containing lower and upper limits
#' @return list of outlisers, values shows their distance from limits
#' @export
HistOutliers2 <- function(div_groups, limits)
{
  upper <- sapply(1:(length(div_groups) - 1), function(w) div_groups[[w]][which(div_groups[[w]] > limits[[w]][[2]])] - limits[[w]][[2]])
  lower <- sapply(1:(length(div_groups) - 1), function(w) div_groups[[w]][which(div_groups[[w]] < limits[[w]][[1]])] - limits[[w]][[1]])
  all_out <- sapply(1:(length(div_groups) - 1), function(w) c(lower[[w]], upper[[w]]))
  return(all_out)
}

#' Does calculation on groups
#'
#' Calculates some statistics defined by functions on list of vectors and returns it as a list of vectors
#' @author Harutyun Khachatryan
#' @param div_groups list of vectors on which statistical calculation should be done
#' @param func_calc name of statistical or other function that can be applied on vector of values
#' @param func_def a function name that should be used in case when func_calc gives NA values, default value is median
#' @return results of calculation for each group
#' @export
CalcOnGroups2 <- function(div_groups, func_calc, func_def = function(w) median(w, na.rm = TRUE))
{
  default_value <- func_def(unlist(div_groups))
  der_vals <- lapply(div_groups, func_calc)
  if (any(is.na(der_vals))) {der_vals[which(is.na(der_vals))] <- default_value}
  return(der_vals)
}
