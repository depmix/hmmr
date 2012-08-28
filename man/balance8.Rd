\name{balance8}

\docType{data}

\alias{balance8}

\title{Balance scale data}

\description{

	Data are from 8 repeated measurements of 25 balance scale items with 1004 
	participants ranging from 6-17 years of age. The data are aggregated over
	5 item types: weight (w), distance (d), conflict weight (cw), conflict 
	distance (cd), and conflict balance (cb). 
	
}

\usage{data(balance8)}

\format{
  A data frame with 168 observations on the following 3 variables.
  \describe{
	\item{\code{rt}}{a numeric vector of response times (log ms)}
	\item{\code{corr}}{a numeric vector of accuracy scores (0/1)}
	\item{\code{Pacc}}{a numeric vector of the pay-off for accuracy}
  }
}

\source{
	
	Brenda Jansen & Ingmar Visser. Manuscript in preparation.
	Corresponding author: i.visser@uva.nl
	
}

\examples{
data(balance8)
## maybe str(balance8) ; plot(balance8) ...
}

\keyword{datasets}
