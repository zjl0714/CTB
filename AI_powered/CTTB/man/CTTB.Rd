\name{CTTB}
\alias{CTTB}
\title{Covariate-tightened trimming bounds}
\description{Estimate aggregated and conditional trimming bounds using covariates.}
\usage{
CTTB(data, seed = NULL, Y, D, S, X = NULL, W = NULL, Pscore = NULL,
     splitting = FALSE, cv_fold = 5, trim_l = 0, trim_u = 1,
     aggBounds = TRUE, IM_cv = TRUE, alpha = 0.05, cBounds = FALSE,
     X_moderator = NULL, LeeBounds = FALSE, cond.mono = TRUE)
}
\arguments{
  \item{data}{A data frame.}
  \item{seed}{Optional random seed used for sample splitting and forest fitting.}
  \item{Y}{Name of the outcome variable.}
  \item{D}{Name of the binary treatment variable.}
  \item{S}{Name of the binary response or selection indicator.}
  \item{X}{Character vector of covariate names.}
  \item{W}{Optional name of a sampling-weight variable. If omitted, equal weights are used.}
  \item{Pscore}{Optional name of a treatment-probability variable. If omitted, propensity scores are estimated from the data.}
  \item{splitting}{Logical; whether to use regression splitting in quantile forests when applicable.}
  \item{cv_fold}{Number of folds used for cross-fitting.}
  \item{trim_l}{Lower trimming threshold applied to estimated probabilities.}
  \item{trim_u}{Upper trimming threshold applied to estimated probabilities.}
  \item{aggBounds}{Logical; whether to estimate aggregated covariate-tightened trimming bounds.}
  \item{IM_cv}{Logical; whether to report Imbens-Manski confidence intervals for the identified set.}
  \item{alpha}{Significance level used for confidence intervals.}
  \item{cBounds}{Logical; whether to estimate conditional bounds.}
  \item{X_moderator}{A data frame or matrix of moderator values at which conditional bounds are evaluated.}
  \item{LeeBounds}{Logical; whether to estimate classic trimming bounds.}
  \item{cond.mono}{Logical; whether to allow conditional monotonicity in the Lee-bounds step.}
}
\details{
The package implements covariate-tightened trimming bounds proposed by Samii,
Wang, and Zhou (2025). The method extends Lee's trimming bounds by using
random-forest estimates of nuisance quantities and cross-fitting to construct
conditional and aggregated identified sets.
}
\value{
A list of class \code{"CTTB"}. Depending on the chosen options, the list can
contain \code{TB}, \code{CTTB}, \code{cB}, \code{parameters}, and \code{call}.
}
\references{
Athey, S., Tibshirani, J., and Wager, S. (2019). Generalized random forests.
\emph{Annals of Statistics}, 47(2), 1148--1178.

Lee, D. S. (2009). Training, wages, and sample selection: Estimating sharp
bounds on treatment effects. \emph{Review of Economic Studies}, 76(3),
1071--1102.

Samii, C., Wang, Y., and Zhou, J. A. (2025). Generalizing trimming bounds for
endogenously missing outcome data using random forests. \emph{Political
Analysis}.
}
\examples{
\dontrun{
data(simData)
res <- CTTB(
  data = simData,
  Y = "Y", D = "D", S = "S",
  X = setdiff(names(simData), c("Y", "D", "S", "Ps")),
  Pscore = "Ps",
  aggBounds = TRUE,
  LeeBounds = TRUE
)
print(res)
}
}
