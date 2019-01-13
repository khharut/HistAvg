\name{periodic}
\alias{periodic}
\title{Periodic dataset}
\description{One month periodic monitor data}
\usage{data(periodic)}
\format{Named list of values, list names are timestamps in miliseconds}
\source{Monitis dashboard}

\keyword{datasets}
\examples{data(periodic)
timestamps<- as.numeric(names(periodic))
values<- as.numeric(periodic)
time_order<- order(timestamps)
timestamps<- timestamps[time_order]
values<- values[time_order]
plot(timestamps, values, type = "l")
}
