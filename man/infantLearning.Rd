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
	A data frame with 532 (19 infants times 28 trials) observations on the
	following 3 variables.  
	
	\describe{ 
	
		item{\code{id}}{infant id}
	
		\item{\code{trial}}{trial number}
		
		\item{center}{the number of looks at the center of the screen.}
		
		\item{away}{the number of looks away from the screen/stimulus.}
				
		\item{acc}{the number of correct anticipatory looks.}

		\item{err}{the number of incorrect anticipatory looks.}

		\item{looks}{the total number of looks during the trial (sum of center,
		away, acc, and err).}
		
		\item{\code{last}}{Whether the last look of infants was a
		'\code{cor}'rect anticipation, an '\code{inc}'orrect antipation or
		'\code{no}' anticipation at all.}

	}
}

\source{
	
	Mandell & Raijmakers (2012), et cetera. 
	
}

\examples{
data(infantLearning)
}

\keyword{datasets}
