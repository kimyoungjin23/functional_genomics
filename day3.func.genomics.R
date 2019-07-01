load("derisi.rda")
source("logratios.derisi.R")
library(GO.db)
keytypes(GO.db)
columns(GO.db)
select(GO.db, keytype="TERM", keys="galactose metabolic process", columns = "GOID")
select(GO.db, keys="GO:0006012", columns=c("TERM","DEFINITION"), keytype="GOID")

library(org.Sc.sgd.db)
keytypes(org.Sc.sgd.db)
columns(org.Sc.sgd.db)

select(org.Sc.sgd.db, keys="GO:0006012", columns="ORF", keytype="GOALL")
GN.table = select(org.Sc.sgd.db, keys="GO:0006012", columns="GENENAME", keytype="GOALL")
all.annotated = keys(org.Sc.sgd.db,keytype="ORF")
length(all.annotated)
select(GO.db, keytype="TERM", keys="tricarboxylic acid cycle", columns = "GOID")
# GO-ID is GO:0006099
sum(pnorm((z7),lower.tail = F)<0.001)
pvals<-which(pnorm((z7),lower.tail = F)<0.001)
tca.GO<-select(org.Sc.sgd.db, keys="GO:0006099", columns="ORF", keytype="GOALL")

orf.is.induced<-universe %in% universe[pvals]

orf.is.TCA<- universe %in% tca.GO$ORF
table(orf.is.induced,orf.is.TCA)

fisher.test(x=orf.is.induced, y=orf.is.TCA) #TCA significant using fisher test

select(GO.db, keytype="TERM", keys="Golgi membrane", columns = "GOID")
gm.GO<-select(org.Sc.sgd.db, keys="GO:0000139", columns="ORF", keytype="GOALL")

orf.is.GM<- universe %in% gm.GO$ORF

table(orf.is.induced,orf.is.GM)
fisher.test(orf.is.induced, orf.is.GM) #NOT significant

N=length(universe)
K=sum(orf.is.TCA)
n=sum(orf.is.induced)
k=sum(universe[orf.is.TCA] %in% universe[orf.is.induced])

phyper(k-1, K, N-K, n, lower.tail = F)# TCA significant based on empirical p value by hypergeometric test
#same results as the fisher test
#the k-1 has to be there in order to account for the hypergeometric distribution that doesn't include k
plot.ecdf(logratios[,7])
plot.ecdf(logratios[,7],add=TRUE,col="blue")
plot.ecdf(logratios[,7][orf.is.TCA],add=TRUE,col="red")
qqplot(logratios[,7],logratios[,7][orf.is.TCA==F])
# they do not overlap because they are enriched differently
######################### t.test for TCA pathway enrichment ################################
logR.not.TCA<-z7[orf.is.TCA==F]
logR.yes.TCA<-z7[orf.is.TCA==T]

stderr1<-sqrt(
  ((length(logR.not.TCA)-1)*var(logR.not.TCA)^2+(length(logR.yes.TCA)-1)*var(logR.yes.TCA)^2)/(length(logR.not.TCA)+length(logR.yes.TCA)-2)
)*sqrt(1/length(logR.not.TCA) + 1/length(logR.yes.TCA))

empirical_tvalue<-(mean(logR.yes.TCA)-mean(logR.not.TCA))/stderr1
 
lograt7box<-boxplot(z7~orf.is.TCA)

myResults<-t.test(z7~ orf.is.TCA)
myResults$statistic
empirical_tvalue

######################### t.test for GM pathway enrichment ################################
plot.ecdf(universe.logratios[,7])
plot.ecdf(universe.logratios[,7][orf.is.GM==F],add=TRUE,col="blue")
plot.ecdf(universe.logratios[,7][orf.is.GM==T],add=TRUE,col="red")
qqplot(logratios[,7],logratios[,7][orf.is.GM==F])

logR.not.GM<-universe.logratios[,7][orf.is.GM==F]
logR.yes.GM<-universe.logratios[,7][orf.is.GM==T]

var_yesno.gm<-sqrt(
  (
    (length(logR.not.GM)-1)*var(logR.not.GM)^2+(length(logR.yes.GM)-1)*var(logR.yes.GM)^2)/(length(logR.not.GM)+length(logR.yes.GM)-2))*sqrt(1/length(logR.not.GM) + 1/length(logR.yes.GM))

empirical_tvalue.gm<(-mean(logR.yes.GM)+mean(logR.not.GM))/var_yesno.gm
empirical_tvalue.gm2<-(-mean(logR.yes.GM)+mean(logR.not.GM))/myResults.gm$stderr

lograt7box.gm<-boxplot(universe.logratios[,7] ~ orf.is.GM)
lograt7box.gm
myResults.gm<-t.test(universe.logratios[,7]~ orf.is.GM)
myResults.gm$statistic
empirical_tvalue.gm
