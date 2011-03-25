
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




hm4 <- depmix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial()), ntimes=rep(8,470), ns=4)
fhm4 <- fit(hm4,emcontrol=em.control(maxit=400))




# do some models ...
library(depmixS4)


# lca models on weight and distance items only
m1 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1),data=d2, family=list(binomial(),binomial()),ns=1)
fm1 <- fit(m1)

m2 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1),data=d2, family=list(binomial(),binomial()),ns=2)
fm2 <- fit(m2)


fm3 <- list()
m3 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1),data=d2, family=list(binomial(),binomial()),ns=3)
for(i in 1:10) {
	fm3[[i]] <- fit(m3)
}


fm3 <- fm3[[1]]

m4 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1),data=d2, family=list(binomial(),binomial()),ns=4)
fm4 <- fit(m4,emcontrol=em.control(maxit=400))
# fm4 <- fit(fm4)

m5 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1),data=d2, family=list(binomial(),binomial()),ns=5)
fm5 <- fit(m5,emcontrol=em.control(maxit=400))

fm5 <- fit(fm5,emcontrol=em.control(maxit=400,random.start=FALSE))

save(fm1,fm2,fm3,fm4,fm5,file="GplusAmodels.Rda")


# lca models on weigth, distance, and conflict weight

# lca models on weight and distance items only
m1 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1),
	data=d2, family=list(binomial(),binomial(),binomial()),ns=1)
fm1 <- fit(m1)

m1 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1),
	data=d2, family=list(binomial(),binomial(),binomial()),ns=2)
fm2 <- fit(m1)

m1 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1),
	data=d2, family=list(binomial(),binomial(),binomial()),ns=3)
fm3 <- fit(m1)

m1 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1),
	data=d2, family=list(binomial(),binomial(),binomial()),ns=4)
fm4 <- fit(m1)

m1 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1),
	data=d2, family=list(binomial(),binomial(),binomial()),ns=5)
fm5 <- fit(m1)

m1 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1),
	data=d2, family=list(binomial(),binomial(),binomial()),ns=6)
fm6 <- fit(m1)




# # use multinomial ipv binomial
# mm2 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1),
# 	data=d2, family=list(multinomial("identity"),multinomial("identity")), ns=2)
# 
# fmm2 <- fit(mm2)



# lca models on all items

ma1 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ns=1)
fma1 <- fit(ma1)


ma2 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ns=2)
fma2 <- fit(ma2,emcontrol=em.control(maxit=400))


ma3 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ns=3)
fma3 <- fit(ma3,emcontrol=em.control(maxit=400))


ma4 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ns=4)
fma4 <- fit(ma4,emcontrol=em.control(maxit=400))


ma5 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ns=5)
fma5 <- fit(ma5,emcontrol=em.control(maxit=400))


ma6 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ns=6)
fma6 <- fit(ma6,emcontrol=em.control(maxit=400))


ma7 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ns=7)
fma7 <- fit(ma7,emcontrol=em.control(maxit=400))


ma8 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ns=8)
fma8 <- fit(ma8,emcontrol=em.control(maxit=400))


ma9 <- mix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ns=9)
fma9 <- fit(ma9,emcontrol=em.control(maxit=400))


save(fma1,fma2,fma3,fma4,fma5,fma6,fma7,fma8,file="lcaAllItems.Rda")


# hmm models on all items

hm1 <- depmix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ns=1, ntimes=rep(8,608))
fhm1 <- fit(hm1,emcontrol=em.control(maxit=400))


hm2 <- depmix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ntimes=rep(8,608), ns=2)
fhm2 <- fit(hm2,emcontrol=em.control(maxit=400))



hm5 <- depmix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ntimes=rep(8,608), ns=5)
fhm5 <- fit(hm5,emcontrol=em.control(maxit=400))


hm6 <- depmix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ntimes=rep(8,608), ns=6)
fhm6 <- fit(hm6,emcontrol=em.control(maxit=400))


hm7 <- depmix(list(cbind(gc,gi)~1,cbind(ac,ai)~1,cbind(cgc,cgi)~1,cbind(cac,cai)~1,cbind(cbc,cbi)~1), 
	data=d2, family=list(binomial(),binomial(),binomial(),binomial(),binomial()), ntimes=rep(8,608), ns=7)
fhm7 <- fit(hm7,emcontrol=em.control(maxit=400))


save(fhm5,fhm6, file="hmmAllItems.Rda")


