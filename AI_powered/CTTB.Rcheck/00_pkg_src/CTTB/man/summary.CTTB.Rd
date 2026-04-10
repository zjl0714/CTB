\name{summary.CTTB}
\alias{summary.CTTB}
\title{Summarize a CTTB object}
\usage{
\method{summary}{CTTB}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{"CTTB"}.}
  \item{...}{Additional arguments, currently unused.}
}
\examples{
\dontrun{
data(simData)
res <- CTTB(simData, Y = "Y", D = "D", S = "S",
            X = setdiff(names(simData), c("Y", "D", "S", "Ps")),
            Pscore = "Ps", aggBounds = TRUE, LeeBounds = TRUE)
summary(res)
}
}
