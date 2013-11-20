\name{infantLearning}

\docType{data}

\alias{infantLearning}

\title{Infant categorization learning data}

\description{

	Data are from 19 infants (11.5-12.5 months) learning to categorize
	objects based on a single feature, either form or color.  Each infant
	performed 28 trials of an occluder paradigm: the stimulus appears,
	disappears behind an occluder, and then reappears at either the left or
	right hand side of the screen.  Learning is determined by anticipatory
	looking to the correct side of the screen.  Importantly, these data also
	contain a measure of habituation: the number of times infants look away
	from the screen on which the stimulus is projected.  Infants tend to
	start paying less attention when stimuli are presented repeatedly, which
	may interfere with determining whether learning has taken place.  The
	data consists of 8 variables describing several aspects of the learning
	process, see details below. 
	
}

\usage{data(infantLearning)}

\format{
  A data frame with 532 observations on the following 8 variables.
  \describe{
    \item{\code{id}}{infant id}
    \item{\code{trial}}{trial number}
    \item{\code{center}}{the number of looks at the center of the screen.}
    \item{\code{away}}{the number of looks away from the screen/stimulus.}
    \item{\code{cor}}{the number of correct anticipatory looks.}
    \item{\code{inc}}{the number of incorrect anticipatory looks.}
    \item{\code{looks}}{the total number of looks during the trial (sum of center,
		away, cor, and inc).}
		
	\item{\code{last}}{a factor with levels \code{cor} \code{inc}
	\code{no}; whether the last look of infants was a '\code{cor}'rect
	anticipation, an '\code{inc}'orrect antipation or '\code{no}'
	anticipation at all.}
		
  }
}

\source{
	
	Mandell & Raijmakers (2012), et cetera. 
	
}

\examples{
data(infantLearning)
}

\keyword{datasets}
