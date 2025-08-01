\name{CTTB}
\alias{CTTB}
\title{Covariate-tightened trimming bounds}
\description{Estimate aggregated and conditional trimming bounds using information from covariates, based on the algorithm in Samii, Wang, and Zhou (2025)}
\usage{CTTB(data, seed = NULL, Y, D, S, X = NULL, W = NULL,
Pscore = NULL, splitting = FALSE, cv_fold = 5,
trim_l = 0, trim_u = 1, aggBounds = TRUE, IM_cv = TRUE,
alpha = 0.05, cBounds = FALSE, X_moderator = NULL, LeeBounds = FALSE,
direction = NULL, cond.mono = TRUE, message = TRUE)
}
\arguments{
\item{data}{a data frame.}

\item{seed}{random seeds used in sampling splitting.}

\item{Y}{name of the outcome variable.}

\item{D}{name of the treatment variable. This variable should be dichotomous (0 and 1) in data.}

\item{S}{name of the variable indicating the response. This variable should be dichotomous (0 and 1) in data.}

\item{X}{names of the covariates.}

\item{W}{name of the variable indicating the sampling weights. If unspecified, units will be equally weighted.}

\item{Pscore}{name of the variable indicating the probability of being treated. If unspecified, the probability will be estimated from data.}

\item{splitting}{a logical flag indicating whether to use regression splitting}

\item{cv_fold}{an integer specifying the number of folds used in cross-validation. The default value is 5.}

\item{trim_l}{the lower bound of }

\item{trim_u}{the upper bound of}

\item{aggBounds}{a logical flag indicating whether to calculate the aggregated bounds. Default is \code{aggBounds = TRUE}.}

\item{IM_cv}{a logical flag indicating whether to use the critical values from Imbens and Manski (2009) for the effect estimate. Default is \code{IM_cv = TRUE}.}

\item{alpha}{the significance level for hypothesis tests and confidence intervals. Default is \code{alpha = 0.05}.}

\item{cBounds}{a logical flag indicating whether to calculate the aggregated bounds. Default is \code{cBounds = TRUE}.}

\item{X_moderator}{covariates values on which the conditional bounds are estimated. If unspecified, the conditional bounds will not be estimated.}

\item{LeeBound}{a logical flag indicating whether to calculate the basic trimming bounds (Lee bounds). Default is \code{LeeBounds = FALSE}.}

\item{direction}{a logical flag indicating whether to calculate the basic trimming bounds (Lee bounds). Default is \code{LeeBounds = FALSE}.}

\item{cond.mono}{a logical flag indicating whether to allow for the assumption of conditional monotonic selection. Default is \code{cond.mono = TRUE}.}
}

\details{
  The \code{CTTB} package implements the covariate-tightened trimming bounds (CTTB) method proposed by Samii, Wang, and Zhou (2025). CTTB generalizes the classic trimming bounds (TB) approach of Lee (2009) by leveraging the generalized random forest algorithm (Athey, Tibshirani, and Wager, 2019) to incorporate information from a potentially high-dimensional set of covariates. The method enables users to construct trimming bounds for the average treatment effect
among always-responders conditional on covariate values, and to aggregate these conditional bounds into a tighter identified set. It is built upon the idea of double machine learning (Chernozhukov et al., 2018) and provides users with honest confidence intervals.
}

\value{
  \item{TB}{estimate of the lower bound in the classic trimming bounds.}

  \item{CTTB}{estimate of the upper bound in the classic trimming bounds.}

  \item{cB}{standard error estimate of the estimated lower bound in the classic trimming bounds.}

  \item{parameters}{a vector of parameters saving.}
}
\author{
  Cyrus Samii (\email{cds2083@nyu.edu}), New York University.

  Ye Wang (\email{yewang@unc.edu}), UNC-Chapel Hill.

  Aaron Junlong Zhou (\email{jlzhou@nyu.edu}), Tecent America
}
\references{
  Athey, Susan, Julie Tibshirani, and Stefan Wager. (2019).
  Generalized Random Forests.
  \emph{The Annals of Statistics, 47}(2): 1148–78.

  Chernozhukov, Victor, Denis Chetverikov, Mert Demirer, Esther Duflo, Christian Hansen, Whitney Newey, and James Robins. (2018).
  Double/Debiased Machine Learning for Treatment and Structural Parameters.
  \emph{The Econometrics Journal, 21}(1): C1–68.

  Lee, David S. (2009).
  Training, Wages, and Sample Selection: Estimating Sharp Bounds on Treatment Effects.
  \emph{The Review of Economic Studies, 76}(3): 1071–1102.

  Samii, Cyrus, Ye Wang, and Junlong Aaron Zhou. (2025).
  Generalizing Trimming Bounds for Endogenously Missing Outcome Data Using Random Forests.
  \emph{Political Analysis, First View}: 1-15.
}

\examples{
library(CTTB)
data(simData)
out <- CTB(data = dat, seed = 1234, Y = "Y", D = "D", S = "S", X = c(names(dat)[c(2:(P+2))]),                       Pscore = "Ps", regression.splitting = FALSE, cv_fold = 5, aggBounds = TRUE,
           cBounds = TRUE, X_moderator = X_moderator, direction = NULL, cond.mono = 0)
print(out)
}


