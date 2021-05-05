\name{balance8models}

\docType{data}

\alias{balance8models}

\title{Models for the balance8 data set}

\description{

	Models with 3-8 states for the balance8 data. 
	
}

\usage{data(balance8models)}
\format{
  A list with fitted models with 3-8 states. 

}

\details{

The dataset is generated with the following code: 

data(balance8)
    \preformatted{
multstart <- function(model, nr=10, initIters=10, verbose=TRUE) {
  llbest <- as.numeric(logLik(model))
  bestmodel <- model
  for(i in 1:nr) {
    fmod <- fit(model, emcontrol=em.control(maxit=initIters)) 
	if(verbose) print(paste(i,": ", logLik(fmod)))
    if(logLik(fmod) > llbest) {
      llbest <- logLik(fmod)
      bestmodel <- fmod
    }
  }
  bestmodel <- fit(bestmodel, emcontrol=em.control(random.start=FALSE))
  return(bestmodel)
}

set.seed(12)

hm3id <- depmix(list(cbind(wc,wi)~1,cbind(dc,di)~1,cbind(cwc,cwi)~1,cbind(cdc,cdi)~1,cbind(cbc,cbi)~1), 
	data=balance8, family=list(binomial("identity"),binomial("identity"),binomial("identity"),binomial("identity"),binomial("identity")), 
	ntimes=rep(8,1004), ns=3,
	respst=rep(0.5,15))
	
fhm3id <- multstart(hm3id)

hm4id <- depmix(list(cbind(wc,wi)~1,cbind(dc,di)~1,cbind(cwc,cwi)~1,cbind(cdc,cdi)~1,cbind(cbc,cbi)~1), 
	data=balance8, family=list(binomial("identity"),binomial("identity"),binomial("identity"),binomial("identity"),binomial("identity")), 
	ntimes=rep(8,1004), ns=4,
	respst=rep(0.5,20))
	
fhm4id <- multstart(hm4id)


hm5id <- depmix(list(cbind(wc,wi)~1,cbind(dc,di)~1,cbind(cwc,cwi)~1,cbind(cdc,cdi)~1,cbind(cbc,cbi)~1), 
	data=balance8, family=list(binomial("identity"),binomial("identity"),binomial("identity"),binomial("identity"),binomial("identity")), 
	ntimes=rep(8,1004), ns=5,
	respst=rep(0.5,25))
	
fhm5id <- multstart(hm5id)

hm6id <- depmix(list(cbind(wc,wi)~1,cbind(dc,di)~1,cbind(cwc,cwi)~1,cbind(cdc,cdi)~1,cbind(cbc,cbi)~1), 
	data=balance8, family=list(binomial("identity"),binomial("identity"),binomial("identity"),binomial("identity"),binomial("identity")), 
	ntimes=rep(8,1004), ns=6,
	respst=rep(0.5,30))
	
fhm6id <- multstart(hm6id)

hm7id <- depmix(list(cbind(wc,wi)~1,cbind(dc,di)~1,cbind(cwc,cwi)~1,cbind(cdc,cdi)~1,cbind(cbc,cbi)~1), 
	data=balance8, family=list(binomial("identity"),binomial("identity"),binomial("identity"),binomial("identity"),binomial("identity")), 
	ntimes=rep(8,1004), ns=7,
	respst=rep(0.5,35))
	
fhm7id <- multstart(hm7id)

set.seed(1)

hm8id <- depmix(list(cbind(wc,wi)~1,cbind(dc,di)~1,cbind(cwc,cwi)~1,cbind(cdc,cdi)~1,cbind(cbc,cbi)~1), 
	data=balance8, family=list(binomial("identity"),binomial("identity"),binomial("identity"),binomial("identity"),binomial("identity")), 
	ntimes=rep(8,1004), ns=8,
	respst=rep(0.5,40))
	
fhm8id <- multstart(hm8id)

balance8models <- list(fhm3id,fhm4id,fhm5id,fhm6id,fhm7id,fhm8id)
}
}


\examples{
data(balance8models)
balance8models[[3]]
}

\keyword{datasets}
