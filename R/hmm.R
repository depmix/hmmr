lca <- function(data, nclasses, fit=TRUE, family=NULL, ...) {
	
	# univariate numeric data by default treated as gaussian
 	if(is.vector(data)) {
		nm <- deparse(substitute(data))
		form <- as.formula(paste(nm,"~1",sep=""))
		if(is.null(family)) family=gaussian()
		simple <- 0
		data <- data.frame(nm=data)
	} else {
		# univariate factor data treated as multinomial
 		if(is.factor(data)) {
			nm <- deparse(substitute(data))
			form <- as.formula(paste(nm,"~1",sep=""))
			if(is.null(family)) family=multinomial("identity")
			simple <- 1
			data <- data.frame(nm=data)
		} else {
			# multivariate data
			form <- list()
			fam <- list()
			nc <- ifelse(!is.null(dim(data)), dim(data)[2],1)
			simple <- numeric(nc)
			for(i in 1:nc) {
				form[[i]] <- as.formula(paste(names(data)[i],"~1",sep=""))
				fam[[i]] <- switch(class(data[,i]),
					"numeric"=gaussian(),
					"factor"=multinomial("identity"),
					stop("Provide family arguments of data other than 'numeric' or 
						'factor's"))
				simple[i] <- switch(class(data[,i]),
					"numeric"=0,
					"factor"=1)
			}
			if(is.null(family)) family=fam
			data <- data
		}
	}
		
	mod <- mix(response=form,data=data,nstates=nclasses,family=family, ...)
	attr(mod,"type") <- "lca"
	res <- fit(mod)
	return(res)
}


setMethod("show","mix",
	function(object) {
		cat("Initial state probabilties model \n")
		print(object@prior)
		cat("\n")
		for(i in 1:object@nstates) {
			cat("Response model(s) for state", i,"\n\n")
			for(j in 1:object@nresp) {
				cat("Response model for response",j,"\n")
				print(object@response[[i]][[j]])
				cat("\n")
			}
			cat("\n")
		}
	}
)

setMethod("summary","mix.fitted",
	function(object,which="all",compact=FALSE) {
		ns <- nstates(object)
		ans=switch(which,
			"all" = 1,
			"response" = 2,
			"prior" = 3,
			stop("Invalid 'which' argument in summary of fitted mix model")
		)
		if(ans==1|ans==3) {
			if(is.null(attr(object,"type"))) { 
				cat("Mixture probabilities model \n")
				print(object@prior)
				cat("\n")
			} else {
				if(attr(object,"type")=="lca") {
					cat("Mixture probabilities \n")
					print(object@prior@parameters$coefficients[1:ns])
				}
			}
		}
		if(ans==1|ans==2) {
			if(!compact) {
				for(i in 1:ns) {
					cat("Response model(s) for state", i,"\n\n")
					for(j in 1:object@nresp) {
						cat("Response model for response",j,"\n")
						print(object@response[[i]][[j]])
						cat("\n")
					}
					cat("\n")
				}
			} else {
				cat("\nResponse parameters \n")
				pars <- list()
				np <- numeric(object@nresp)
				for(j in 1:object@nresp) {
					np[j] <- npar(object@response[[1]][[j]])
					pars[[j]] <- matrix(,nr=ns,nc=np[j])
				}
				allpars <- matrix(,nr=ns,nc=0)
				for(j in 1:object@nresp) {
					for(i in 1:ns) {
						pars[[j]][i,]=getpars(object@response[[i]][[j]])
					}
					allpars <- cbind(allpars,pars[[j]])
				}
				rownames(allpars) <- paste("St",1:ns,sep="")
				print(allpars)
			}
		}
	}	
)


hmm <- function(data, nstates, fit=TRUE, ntimes=NULL, family=NULL) {
	
	# univariate numeric data by default treated as gaussian
 	if(is.vector(data)) {
		nm <- deparse(substitute(data))
		form <- as.formula(paste(nm,"~1",sep=""))
		if(is.null(family)) family=gaussian()
		if(is.null(ntimes)) ntimes <- length(data)
		simple <- 0
		data <- NULL
	} else {
		# univariate factor data treated as multinomial
 		if(is.factor(data)) {
			nm <- deparse(substitute(data))
			form <- as.formula(paste(nm,"~1",sep=""))
			if(is.null(ntimes)) ntimes <- length(data)
			if(is.null(family)) family=multinomial("identity")
			simple <- 1
			data <- NULL
		} else {
			# multivariate data
			form <- list()
			fam <- list()
			nc <- ifelse(!is.null(dim(data)), dim(data)[2],1)
			simple <- numeric(nc)
			for(i in 1:nc) {
				form[[i]] <- as.formula(paste(names(data)[i],"~1",sep=""))
				fam[[i]] <- switch(class(data[,i]),
					"numeric"=gaussian(),
					"factor"=multinomial("identity"),
					stop("Provide family arguments of data other than 'numeric' or 
						'factor's"))
				simple[i] <- switch(class(data[,i]),
					"numeric"=0,
					"factor"=1)
			}
			if(is.null(family)) family=fam
			if(is.null(ntimes)) ntimes=nrow(data)
			data <- data
		}
	}
		
	mod <- depmix(response=form,data=data,nstates=nstates,ntimes=ntimes,family=family)
	attr(mod,"type") <- "hmm"
	res <- fit(mod)
	return(res)
}


# 
# PRINT method
# 

setMethod("show","depmix",
	function(object) {
		ns <- nstates(object)
		if(is.null(attr(object,"type"))) { 
			cat("Initial state probabilties model \n")
			print(object@prior)
			cat("\n")
			for(i in 1:object@nstates) {
				cat("Transition model for state (component)", i,"\n")
				print(object@transition[[i]])
				cat("\n")
			}
			cat("\n")
		} else {
			if(attr(object,"type")=="hmm") { 
				cat("Initial state probabilties\n")
				print(object@prior@parameters$coefficients[1:ns])
				cat("\nTransition matrix \n")
				pars <- getpars(object)
				trm <- matrix(pars[(ns+1):(ns^2+ns)],ns,ns,byr=T)
				rownames(trm) <- paste("fromS",1:ns,sep="")
				colnames(trm) <- paste("toS",1:ns,sep="")
				print(trm)
				cat("\n")
			} 
		}
		for(i in 1:ns) {
			cat("Response model(s) for state", i,"\n\n")
			for(j in 1:object@nresp) {
				cat("Response model for response",j,"\n")
				print(object@response[[i]][[j]])
				cat("\n")
			}
			cat("\n")
		}
	}
)

setMethod("summary","depmix.fitted",
	function(object,which="all", compact=FALSE) {
		ns <- object@nstates
		ans=switch(which,
			"all" = 1,
			"response" = 2,
			"prior" = 3,
			"transition" = 4,
			stop("Invalid 'which' argument in summary of fitted depmix model")
		)
		if(ans==1|ans==3) {
			if(is.null(attr(object,"type"))) { 
				cat("Initial state probabilties model \n")
				print(object@prior)
				cat("\n")
			} else {
				if(attr(object,"type")=="hmm") {
					cat("Initial state probabilties\n")
					print(object@prior@parameters$coefficients[1:ns])
				}
			}
		}
		if(ans==1|ans==4) {
			if(is.null(attr(object,"type"))) { 
				for(i in 1:ns) {
					cat("Transition model for state (component)", i,"\n")
					print(object@transition[[i]])
					cat("\n")
				}
				cat("\n")
			} else {
				if(attr(object,"type")=="hmm") {
					cat("\nTransition matrix \n")
					pars <- getpars(object)
					trm <- matrix(pars[(ns+1):(ns^2+ns)],ns,ns,byr=T)
					rownames(trm) <- paste("fromS",1:ns,sep="")
					colnames(trm) <- paste("toS",1:ns,sep="")
					print(trm)
					cat("\n")
				}
			}
		}
		if(ans==1|ans==2) {
			if(!compact) {
				for(i in 1:ns) {
					cat("Response model(s) for state", i,"\n\n")
					for(j in 1:object@nresp) {
						cat("Response model for response",j,"\n")
						print(object@response[[i]][[j]])
						cat("\n")
					}
					cat("\n")
				}
			} else {
				cat("Response parameters \n")
				pars <- list()
				np <- numeric(object@nresp)
				for(j in 1:object@nresp) {
					np[j] <- npar(object@response[[1]][[j]])
					pars[[j]] <- matrix(,nr=ns,nc=np[j])
				}
				allpars <- matrix(,nr=ns,nc=0)
				for(j in 1:object@nresp) {
					for(i in 1:ns) {
						pars[[j]][i,]=getpars(object@response[[i]][[j]])
					}
					allpars <- cbind(allpars,pars[[j]])
				}
				rownames(allpars) <- paste("St",1:ns,sep="")
				print(allpars)
			}
		}
	}
)

