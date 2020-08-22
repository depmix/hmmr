\name{MNAR_simulation_results}

\alias{MNAR_simulation_results}

\docType{data}

\title{
  Missing not at random (MNAR) simulation results.
}

\description{
  Results of a simulation study on the effect of missing data on estimates of parameters of a hidden Markov model. In the simulation, 1000 datasets were simulated, each consisting of 100 time series of length 50. The model generating the data was a hidden Markov model with 3 states. Missing values were generated for 5\% of cases in state 1, 25\% of cases in state 2, and 50\% of cases in state 3 (i.e. data can be considered missing not at random, as missingness is state-dependent).
}

\usage{data("NMAR_simulation_results")}
\format{
  A \code{data.frame} with the following variables:
  \describe{
    \item{\code{parameter}}{Parameter name}
    \item{\code{true}}{True value of parameter}
    \item{\code{MAR_estimates_mean}}{Average of parameter estimates for a model that assumes MAR.}
    \item{\code{MAR_estimates_sd}}{Standard deviation of parameter estimates for a model that assumes MAR.}
    \item{\code{MAR_estimates_MAE}}{Mean Absolute Error of parameter estimates for a model that assumes MAR.}
    \item{\code{MNAR_estimates_mean}}{Average of parameter estimates for a model that assumes MNAR.}
    \item{\code{MNAR_estimates_sd}}{Standard deviation of parameter estimates for a model that assumes MNAR.}
    \item{\code{MNAR_estimates_MAE}}{Mean Absolute Error of parameter estimates for a model that assumes MNAR.}
    \item{\code{percMAE}}{Relative MAE of the MNAR model over the MAR model (e.g. MAE_MNAR/MAE_MAR)}
  }
}

\examples{
  data(MNAR_simulation_results)
}
\keyword{datasets}
