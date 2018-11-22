\name{dccs}

\docType{data}

\alias{dccs}

\title{Dimensional Change Card Sort Task Data}

\description{
	
	Data from 3-5 year old children performing on the dimensional 
	change card sort task published in Van Bers et al (2011). 
	
}

\usage{
data(dccs)
}

\format{

A data.frame consisting of the following variable:
  \describe{
	\item{\code{pp}}{participant number}
  	\item{\code{ageM}}{(numeric) age in months}
  	\item{\code{ageY}}{(numeric) age in years}
  	\item{\code{ageGr}}{(numeric) age group}
  	\item{\code{sex}}{(numeric) age group}
  	\item{\code{nTrPre}}{(numeric) number of trials completed pre-switch}
  	\item{\code{nTrPre}}{(numeric) number of trials correct in pre-switch trials}
  	\item{\code{t1Post}}{(numeric) indicator variable for post switch trial 
	indicating correct (1) or incorrect (0)}
  	\item{\code{t2Post}}{(numeric) indicator variable for post switch trial 
	indicating correct (1) or incorrect (0)}
  	\item{\code{t3Post}}{(numeric) indicator variable for post switch trial 
	indicating correct (1) or incorrect (0)}
  	\item{\code{t4Post}}{(numeric) indicator variable for post switch trial 
	indicating correct (1) or incorrect (0)}
  	\item{\code{t5Post}}{(numeric) indicator variable for post switch trial 
	indicating correct (1) or incorrect (0)}
  	\item{\code{t6Post}}{(numeric) indicator variable for post switch trial 
	indicating correct (1) or incorrect (0)}
  \item{\code{nCorPost}}{(numeric) number of trials correct post-switch}
  \item{\code{passPost}}{(numeric) indicator variable whether post-switch test was passed}
  }
  
}

\source{

	Bianca M.C.W. van Bers, Ingmar Visser, Tessa J.P. van Schijndel, 
	Dorothy J. Mandell and Maartje E.J. Raijmakers (2011). The dynamics 
	of development on the Dimensional Change Card Sorting task. 
	\emph{Developmental Science, vol 14(5)}, 960-971.
	
}

\keyword{datasets}
