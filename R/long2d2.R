
setwd("~/Documents/Dropbox/hmmr/datasets")

dt <- read.spss("long.sav",to.data.frame=TRUE)

attach(dt)

sel <- c("NUMMER","GESLACHT","GEBOORTE", "SCHOOL",   "GROEP",    
	"GROEPSL", "DAYS",     "YEARS",    
	"AANWTM1", "AANWTM2",  "AANWTM3",   "AANWTM4" ,
	"AANWTM5",  "AANWTM6",  "AANWTM7", "AANWTM8", 
	"REGELS1",  "REGELS2", "REGELS3",  "REGELS4",  "REGELS5",  "REGELS6",  
	"REGELS7",  "REGELS8", 
	"SCLLEVEL")


# select relevant columns
d1 <- dt[,sel]

d2 <- reshape(d1,varying=list(c("AANWTM1", "AANWTM2",  "AANWTM3",   "AANWTM4" ,
	"AANWTM5",  "AANWTM6",  "AANWTM7", "AANWTM8"),
	c("REGELS1",  "REGELS2", "REGELS3",  "REGELS4",  "REGELS5",  "REGELS6",  
	"REGELS7",  "REGELS8")),sep="",direction="long")

# compute sum scores of different item types

afn <- paste("AFN",1:8,sep="")

# weight items
afng <- paste(afn,c("_G"),sep="")
afng <- sapply(afng,paste,1:5,sep="")

# afname 1
gc <- apply(dt[,afng[,1]],1,sum,na.rm=TRUE)
gi <- apply(1-dt[,afng[,1]],1,sum,na.rm=TRUE)

# afname 2-8
for(i in 2:8) {
	gc <- c(gc,apply(dt[,afng[,i]],1,sum,na.rm=TRUE))
	gi <- c(gi,apply(1-dt[,afng[,i]],1,sum,na.rm=TRUE))
}

gs <- gc + gi

d2$"gc" <- gc
d2$"gi" <- gi
d2$"gs" <- gs



# distance items
afna <- paste(afn,c("_A"),sep="")
afna <- sapply(afna,paste,1:5,sep="")

# afname 1
ac <- apply(dt[,afna[,1]],1,sum,na.rm=TRUE)
ai <- apply(1-dt[,afna[,1]],1,sum,na.rm=TRUE)

# afname 2-8
for(i in 2:8) {
	ac <- c(ac,apply(dt[,afna[,i]],1,sum,na.rm=TRUE))
	ai <- c(ai,apply(1-dt[,afna[,i]],1,sum,na.rm=TRUE))
}

as <- ac + ai

d2$"ac" <- ac
d2$"ai" <- ai
d2$"as" <- as



# conflict weight items
afncg <- paste(afn,c("_CG"),sep="")
afncg <- sapply(afncg,paste,1:5,sep="")

# afname 1
cgc <- apply(dt[,afncg[,1]],1,sum,na.rm=TRUE)
cgi <- apply(1-dt[,afncg[,1]],1,sum,na.rm=TRUE)

# afname 2-8
for(i in 2:8) {
	cgc <- c(cgc,apply(dt[,afncg[,i]],1,sum,na.rm=TRUE))
	cgi <- c(cgi,apply(1-dt[,afncg[,i]],1,sum,na.rm=TRUE))
}

cgs <- cgc + cgi

d2$"cgc" <- cgc
d2$"cgi" <- cgi
d2$"cgs" <- cgs



# conflict distance
afnca <- paste(afn,c("_CA"),sep="")
afnca <- sapply(afnca,paste,1:5,sep="")

# afname 1
cac <- apply(dt[,afnca[,1]],1,sum,na.rm=TRUE)
cai <- apply(1-dt[,afnca[,1]],1,sum,na.rm=TRUE)

# afname 2-8
for(i in 2:8) {
	cac <- c(cac,apply(dt[,afnca[,i]],1,sum,na.rm=TRUE))
	cai <- c(cai,apply(1-dt[,afnca[,i]],1,sum,na.rm=TRUE))
}

cas <- cac + cai

d2$"cac" <- cac
d2$"cai" <- cai
d2$"cas" <- cas



# conflict balance
afncb <- paste(afn,c("_CB"),sep="")
afncb <- sapply(afncb,paste,1:5,sep="")

# afname 1
cbc <- apply(dt[,afncb[,1]],1,sum,na.rm=TRUE)
cbi <- apply(1-dt[,afncb[,1]],1,sum,na.rm=TRUE)

# afname 2-8
for(i in 2:8) {
	cbc <- c(cbc,apply(dt[,afncb[,i]],1,sum,na.rm=TRUE))
	cbi <- c(cbi,apply(1-dt[,afncb[,i]],1,sum,na.rm=TRUE))
}

cbs <- cbc + cbi

d2$"cbc" <- cbc
d2$"cbi" <- cbi
d2$"cbs" <- cbs


# score op oefenitems

afnoef <- paste(afn,c("_"),sep="")
afnoef <- sapply(afnoef,paste,1:3,sep="")
afnoef <- matrix(paste(afnoef,c("C"),sep=""),3,8)

# afname 1
oef <- apply(dt[,afnoef[,1]],1,sum,na.rm=TRUE)

# afname 2-8
for(i in 2:8) {
	oef <- c(oef,apply(dt[,afnoef[,i]],1,sum,na.rm=TRUE))
}

d2$"oef" <- oef


save(d2, file="d2.Rda")




