\name{WPT}

\docType{data}

\alias{WPT}

\title{Weather Prediction Task Data}

\description{

	This data set is a bivariate series of response times and accuracy
	scores of a single participant switching between slow/accurate
	responding and fast guessing on a lexical decision task.  The slow and
	accurate responding, and the fast guessing can be modelled using two
	states, with a switching regime between them.  The dataset further
	contains a third variable called Pacc, representing the relative
	pay-off for accurate responding, which is on a scale of zero to one.
	The value of Pacc was varied during the experiment to induce the
	switching.  This data set is a from participant A in experiment 
	1a from Dutilh et al (2011). The data here is the first series of
	168 trials. The \code{speed} data set in the \code{depmixS4} package 
	has two more series of 134 and 137 trials respectively. 
	
}

\usage{data(WPT)}

\format{
  A data frame with 9600 observations on the following variables.
  \describe{
  \item{\code{id}}{a factor with participant IDs}
  \item{\code{group}}{a factor with group IDs (Parksinson's patient or control)}
  \item{\code{med}}{a factor indicating, for the PD patients, whether they were on dopaminergic medicine or not}
  \item{\code{occ}}{a numeric vector with testing occassions}
  \item{\code{trial}}{a numeric vector with trial numbers}
  \item{\code{c1}}{a numeric (binary) vector indicating whether the first cue was present (1) or not (0)}
  \item{\code{c2}}{a numeric (binary) vector indicating whether the second cue was present (1) or not (0)}
  \item{\code{c3}}{a numeric (binary) vector indicating whether the third cue was present (1) or not (0)}
  \item{\code{c4}}{a numeric (binary) vector indicating whether the fourth cue was present (1) or not (0)}
  \item{\code{y}}{a factor with the actual outcome (Rainy or Fine)}
  \item{\code{r}}{a factor with participants' prediction of the outcome}
  }
}

\source{

Speekenbrink, M., Lagnado, D. A., Wilkinson, L., Jahanshahi, M., & Shanks, D. R. (2010). Models of probabilistic category learning in Parkinson's disease: Strategy use and the effects of L-dopa. \emph{Journal of Mathematical Psychology}, \emph{54}, 123-136.
	
Corresponding author: m.speekenbrink@ucl.ac.uk
	
}

\examples{
data(WPT)
## maybe str(WPT) ; plot(WPT) ...
}

\keyword{datasets}
