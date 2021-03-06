\name{speed_boot_par}
\alias{speed_boot_par}
\docType{data}
\title{
speed boot par
}
\description{
Example of a parametric bootstrap for paraneter inference
}
\usage{data("speed_boot_par")}
\format{
  A \code{boot} object
}
\details{
The bootstrap sample was generated by the following code:

\preformatted{require(depmixS4)
require(hmmr)
require(boot)

data(speed1)
set.seed(5)
spmix2 <- mix(RT~1, data=speed1, nstates=2)
fspmix2 <- fit(spmix2,verbose=FALSE)

# define a function to produce a bootstrap sample 
speed.rg <- function(data,mle) {
  simulate(data)
}
# define what to do with a sample (i.e. estimate parameters)
speed.fun <- function(data) { 
  getpars(fit(data,verbose=FALSE,emcontrol=
              em.control(random.start=FALSE)))
}
# produce 1000 bootstrap samples (may take some time!)
speed_boot_par <- boot(fspmix2,speed.fun,R=1000,sim="parametric",
                   ran.gen = speed.rg)
}
}
\examples{
data(speed_boot_par)
# confidence intervals
confint <- apply(speed_boot_par$t,2,quantile,probs=c(.025,.975))
colnames(confint) <- c("p1","p2","m1","sd1","m2","sd2")
confint
}
\keyword{datasets}
