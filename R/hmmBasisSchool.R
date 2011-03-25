setwd("~/Documents/Dropbox/hmmr/datasets")

load(file="d2.Rda")

# reorder to arrange data per subject
d2 <- d2[order(d2$"NUMMER"),]

# select complete data
sel <- d2$NUMMER[which(d2$"AANWTM1"=="Selected"&d2$time==8)]
d2 <- d2[d2$NUMMER%in%sel,]

sel2 <- "basis"

d2 <- d2[d2$SCLLEVEL%in%sel2,]

pps <- unique(d2$NUMMER)[which(by(d2$oef,list(d2$NUMMER),sum)>22)]

d2 <- d2[d2$SCLLEVEL%in%pps,]


# onderstaand geeft een mooi model met R0, R1, R2, R3/ADD
# meer startwaarden nodig, beste likelihood tot nu toe: -7788.696
# (voor het 4-state model)

hm4 <- depmix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial()), ntimes=rep(8,470), ns=4)
fhm4 <- fit(hm4,emcontrol=em.control(maxit=400))


# beste ll tot nu toe: -7364
hm5 <- depmix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial()), 
	ntimes=rep(8,470), ns=5)
fhm5 <- fit(hm5,emcontrol=em.control(maxit=400))



hm6 <- depmix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial()), 
	ntimes=rep(8,470), ns=6)
fhm6 <- fit(hm6,emcontrol=em.control(maxit=400))

hm7 <- depmix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial()), 
	ntimes=rep(8,470), ns=7)
fhm7 <- fit(hm7,emcontrol=em.control(maxit=400))
