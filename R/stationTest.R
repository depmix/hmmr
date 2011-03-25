
library(depmixS4)

setwd("~/Documents/projects/depmixProject/codesvn")

setwd("/Users/ivisser1/Documents/Dropbox/hmmfunction")

options(digits=4)


data(speed)
speed$corr <- as.factor(speed$corr)

data(balance)
balance <- balance[order(balance$agedays),]
bal <- balance[,4:6]
bal$t1 <- as.factor(bal$t1)
bal$t2 <- as.factor(bal$t2)
bal$t3 <- as.factor(bal$t3)

source("hmm.R")

# four binary items on the balance scale task
mod4 <- mix(list(t1~1,t2~1,t3~1), data=balance, nstates=2,
	family=list(multinomial("identity"),multinomial("identity"),multinomial("identity")))

fm4 <- fit(mod4)

summary(fm4)

m3 <- lca(bal,nc=2)

m3 <- lca(bal[,1:2],nc=2)



summary(m3,com=T)

h1 <- hmm(speed[,1],ns=2)

h1

summary(h1,com=T)

h2 <- hmm(speed[,1],ns=2)

summary(h2,com=T)


h3 <- hmm(speed$rt,ns=2)
summary(h3,com=T)


h4 <- hmm(speed[,1:2],ns=2)
summary(h4,com=T)


m2 <- lca(speed[,1:2],nc=2)






