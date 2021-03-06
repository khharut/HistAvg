% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HistAvg2.R
\name{HistAvg2}
\alias{HistAvg2}
\title{Dynamic threshold calculator}
\usage{
HistAvg2(json.str, figure.out = FALSE, legacy.plot = FALSE, dmin = NA,
  out.path = " ", sens = "low", test.mode = FALSE)
}
\arguments{
\item{json.str}{JSON input can be both string or file path}

\item{figure.out}{logical TRUE/FALSE indicating if after analysis a plot of analysis will be created}

\item{legacy.plot}{logical TRUE/FALSE indicating if plot will be new style i.e. last three days with thresholds
or old one i.e. 2 plots: 1) raw data plot and 2) all raw data points, last day line and thresholds
(default value is FALSE)}

\item{dmin}{fixed step (minutes) between two neighbours in derived dynamic threshold. Default value is 10, but
in case when it is impossible to calculate with that step, then it autocorrected}

\item{out.path}{plot output file folder path, should be accessible for write. Default value is " ", in that
case plot file will be generated in working directory of R and plot file path anyway is written in output information
generated by this method}

\item{sens}{shows confidence interval has three values "low", "medium" and "high" (default is "low")}

\item{test.mode}{logical TRUE/FALSE showing if test mode is enabled or not, default is FALSE}
}
\value{
JSON string containing
a) field "data" containing [timestamp, lower, upper]... upper and lower thresholds along with timestamps
b) field "quality" containing estimated threshold quality value from 0 to 100
c) field "periodicty" containing periodicty reliability
d) field "image" containing full path of created image when figure.out input paramater is TRUE
e) field "errMsg" containing error messages of any error occured during calculation
f) field "warnMsg" containing warning messages of any non-critical error occured during calculation
}
\description{
Finds dynamic threshold for data in JSON format
}
\examples{
library("HistAvg2")
data(periodic)
jsondata<- RJSONIO::toJSON(periodic)
HistAvg2(jsondata)
}
\author{
Harutyun Khachatryan, Tigran Khachikyan, Anoush Ghambaryan
}
