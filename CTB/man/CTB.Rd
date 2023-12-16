\name{CTB}
\alias{CTB}
\title{Covariate-tightened trimming bounds}
\description{Estimate aggregated and conditional trimming bounds using information from covariates, based on the algorithm in Samii, Wang, and Zhou (2023)}
\usage{CTB(data, seed = NULL, Y, D, S, X = NULL, W = NULL, Pscore = NULL,
           regression.splitting = FALSE, cv_fold = 5, trim_l = 0, trim_u = 1,
           aggBounds = TRUE, cBounds = FALSE, X_moderator = NULL,
           direction = NULL, cond.mono = TRUE)
}
\arguments{
\item{data}{a data frame.}

\item{seed}{random seeds used in sampling splitting.}

\item{Y}{the outcome indicator.}

\item{D}{the treatment indicator. It should be a binary variable (0 and 1)}

\item{S}{the response indicator. It should be a binary variable (0 and 1)}

\item{X}{covariates.}

\item{W}{the indicator for sampling weights. If unspecified, units will be equally weighted.}

\item{Pscore}{the indicator for probability of being treated. If unspecified, the probability will be estimated from data.}

\item{regression.splitting}{a logical flag indicating whether to use regression splitting}

\item{cv_fold}{an integer specifying the number of folds used in cross-validation. The default value is 5.}

\item{trim_l}{}

\item{trim_u}{the upper bound of}

\item{aggBounds}{a logical flag indicating whether to calculate the aggregated bounds. Default is \code{aggBounds = TRUE}.}

\item{cBounds}{a logical flag indicating whether to calculate the aggregated bounds. Default is \code{cBounds = TRUE}.}

\item{X_moderator}{covariates values on which the conditional bounds are estimated. If unspecified, the conditional bounds will not be estimated.}

\item{X_moderator}{a logical flag indicating whether to calculate the aggregated bounds. Default is \code{aggBounds = TRUE}.}

\item{cond.mono}{a logical flag indicating whether to allow for the assumption of conditional monotonic selection. Default is \code{cond.mono = TRUE}.}
}

\details{
  \code{CTB} estimates aggregated and conditional trimming bounds using information from covariates.
}

\value{
  \item{tau_l_est_lee}{estimate of the lower bound in the classic trimming bounds.}

  \item{tau_u_est_lee}{estimate of the upper bound in the classic trimming bounds.}

  \item{se_tau_l_est_lee}{standard error estimate of the estimated lower bound in the classic trimming bounds.}

  \item{se_tau_u_est_lee}{standard error estimate of the estimated upper bound in the classic trimming bounds.}

  \item{tau_l_est_avg}{estimate of the lower bound in the aggregated covariate-tightened trimming bounds.}

  \item{tau_u_est_avg}{estimate of the upper bound in the aggregated covariate-tightened trimming bounds.}

  \item{se_tau_l_est_avg}{standard error estimate of the estimated lower bound in the aggregated covariate-tightened trimming bounds.}

  \item{se_tau_u_est_avg}{standard error estimate of the estimated upper bound in the aggregated covariate-tightened trimming bounds.}

  \item{tau_l_est_cond}{estimates of the lower bounds in the conditional covariate-tightened trimming bounds.}

  \item{tau_u_est_cond}{estimates of the upper bounds in the conditional covariate-tightened trimming bounds.}

  \item{se_tau_l_est_cond}{pointwise standard error estimates of the estimated lower bounds in the conditional covariate-tightened trimming bounds.}

  \item{se_tau_u_est_cond}{pointwise standard error estimates of the estimated upper bounds in the conditional covariate-tightened trimming bounds.}
}
\author{
  Cyrus Samii; Ye Wang; Junlong Aaron Zhou
}
\references{
  Lee, David S. 2003. "Training, wages, and sample selection: Estimating sharp bounds on treatment effects." \emph{Review of Economic Studies} 76.3 (2009): 1071-1102.

  Samii, Cyrus, Ye Wang, and Junlong Aaron Zhou. 2023. "Generalizing Trimming Bounds for Endogenously
Missing Outcome Data Using Random Forests." Working paper. Available at: \url{https://arxiv.org/abs/2309.08985}.
}

\examples{
library(CTB)
data(simData)
out <- CTB(data = dat, seed = 1234, Y = "Y", D = "D", S = "S", X = c(names(dat)[c(2:(P+2))]),                       Pscore = "Ps", regression.splitting = FALSE, cv_fold = 5, aggBounds = TRUE,
           cBounds = TRUE, X_moderator = X_moderator, direction = NULL, cond.mono = 0)
print(out)
}


