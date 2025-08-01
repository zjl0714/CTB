\name{CTTB-package}
\alias{CTTB-package}
\title{Covariate-tightened trimming bounds}
\description{Estimate aggregated and conditional trimming bounds for treatment effects using information from covariates, based on the algorithm in Samii, Wang, and Zhou (2025)}

\details{The \code{CTTB} package implements the covariate-tightened trimming bounds (CTTB) method proposed by Samii, Wang, and Zhou (2025). CTTB generalizes the classic trimming bounds (TB) approach of Lee (2009) by leveraging the generalized random forest algorithm (Athey, Tibshirani, and Wager, 2019) to incorporate information from a potentially high-dimensional set of covariates. The method enables users to construct trimming bounds for the average treatment effect
among always-responders conditional on covariate values, and to aggregate these conditional bounds into a tighter identified set. It is built upon the idea of double machine learning (Chernozhukov et al., 2018) and provides users with honest confidence intervals.
}
\author{
  Cyrus Samii (\email{cds2083@nyu.edu}), New York University
  Ye Wang (\email{yewang@unc.edu}), UNC-Chapel Hill
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

