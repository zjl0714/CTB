\name{print.CTTB}
\alias{print.CTTB}
\title{Print a CTTB object}
\usage{
\method{print}{CTTB}(x, ...)
}
\arguments{
  \item{x}{An object of class \code{"CTTB"}.}
  \item{...}{Additional arguments, currently unused.}
}
\examples{
\dontrun{
data(simData)
res <- CTTB(simData, Y = "Y", D = "D", S = "S",
            X = setdiff(names(simData), c("Y", "D", "S", "Ps")),
            Pscore = "Ps", aggBounds = TRUE, LeeBounds = TRUE)
print(res)
}
}
