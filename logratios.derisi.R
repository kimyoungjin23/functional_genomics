load("derisi.rda")
library("org.Sc.sgd.db")
logratios = log2(derisi[,3:9]-derisi[,10:16])-log2(derisi[,17:23]-derisi[,24:30])
rownames(logratios)=derisi$ORF

all.annotated = keys(org.Sc.sgd.db,keytype="ORF")
universe = intersect(all.annotated,rownames(logratios))
logratios = logratios[universe,]

m0 = mean(logratios$R1)
sd0 = sd(logratios$R1)
z1 = (logratios$R1-m0)/sd0
z7 = (logratios$R7-m0)/sd0
p1 = pnorm(z1,lower.tail = F)
p7 = pnorm(z7,lower.tail = F)
universe[which.min(p7)]
sum(p7<0.001) # number of "discoveries"
length(p7)*0.001 #how many p values would be less than 0.001 
