\name{plot.CTTB}
\alias{plot.CTTB}
\title{Plot CTTB results}
\usage{
\method{plot}{CTTB}(x, y = NULL, type = "aggregated", moderator_plot = 1, ...)
}
\arguments{
  \item{x}{An object of class \code{"CTTB"}.}
  \item{y}{Unused.}
  \item{type}{Either \code{"aggregated"} or \code{"conditional"}.}
  \item{moderator_plot}{Index of the moderator to plot when \code{type = "conditional"}.}
  \item{...}{Additional arguments, currently unused.}
}
\value{A \code{ggplot2} object.}
\examples{
\dontrun{
data(simData)
res <- CTTB(simData, Y = "Y", D = "D", S = "S",
            X = setdiff(names(simData), c("Y", "D", "S", "Ps")),
            Pscore = "Ps", aggBounds = TRUE, LeeBounds = TRUE)
plot(res)
}
}
