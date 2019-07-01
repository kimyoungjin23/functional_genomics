ru<-runif(1000,0,1)
hist(ru)
as.matrix(ru)
rm<-matrix(ru, ncol=100)
rm.mean=apply(rm, 1, mean)
hist(rm.mean)
rm.qt=apply(rm, 1, quantile,0.95)
hist(rm.qt)

# https://www.yeastgenome.org/run_seqtools?format=fasta&type=coding&genes=YAL001C&strains=S288C
urlpieces<-vector()
urlpieces[1]<-c("https://www.yeastgenome.org/run_seqtools?format=fasta&type=coding&genes=")
urlpieces[2]<-c("&strains=S288C")
yeastFasta = function (orf) {
  print(orf)
  url=paste0(urlpieces[1],orf,urlpieces[2])
  print(paste0("url=",url))
  readDNAStringSet(url)
  }

BiocManager::install("Biostrings")
library("Biostrings")

orfs=c("YHR023W", "YAL001C")
y=sapply(orfs,yeastFasta)

source("readFasta.R")
BiocManager::install("AnnotationHub")
library("AnnotationHub")
ah=AnnotationHub()

dm=display(ah)

sub_ah = query(ah, c("OrgDb", "cerevisiae"))
sub_ah
unique(sub_ah$species)
orgdb <- query(sub_ah, "OrgDb")[[1]]
orgdb
head(keys(orgdb, keytype="ORF"))
head(keys(orgdb, keytype="GO"))
orfids = c("YAL001C", "YAL002W", "YAL003W", "YAL004W", "YAL005C")
select(orgdb, keys = orfids, columns = "GENENAME", keytype="ORF")

sub_ah2<-query(ah, "TxDb")
txdb<-query(sub_ah2,"hg38")


BiocManager::install("GenomicFeatures")
library("GenomicFeatures")
txdb = ah[["AH52260"]]
tx<-transcripts(txdb)

ah_roadmap<-query(ah, c("RoadMap", "narrowpeak"))
length(ah_roadmap)
head(ah_roadmap)
query(ah_roadmap, "DNAse")

load("derisi.rda")
BiocManager::install("GO.db")
BiocManager::install("org.Sc.sgd.db")
library("GO.db")
library("org.Sc.sgd.db")
source("First_day.R")

load("derisi.rda")
ACR1<-which(derisi$Name == "ACR1")
dim(derisi)
derisi$ORF[which(derisi$Name == "ACR1")]
intensities.R<-matrix(data=NA,nrow=nrow(derisi),ncol=7)
intensities.R[,1]<-derisi$R1-derisi$R1.Bkg
intensities.R[,2]<-derisi$R2-derisi$R2.Bkg
intensities.R[,3]<-derisi$R3-derisi$R3.Bkg
intensities.R[,4]<-derisi$R4-derisi$R4.Bkg
intensities.R[,5]<-derisi$R5-derisi$R5.Bkg
intensities.R[,6]<-derisi$R6-derisi$R6.Bkg
intensities.R[,7]<-derisi$R7-derisi$R7.Bkg

intensities.G<-matrix(data=NA,nrow=nrow(derisi),ncol=7)
intensities.G[,1]<-derisi$G1-derisi$G1.Bkg
intensities.G[,2]<-derisi$G2-derisi$G2.Bkg
intensities.G[,3]<-derisi$G3-derisi$G3.Bkg
intensities.G[,4]<-derisi$G4-derisi$G4.Bkg
intensities.G[,5]<-derisi$G5-derisi$G5.Bkg
intensities.G[,6]<-derisi$G6-derisi$G6.Bkg
intensities.G[,7]<-derisi$G7-derisi$G7.Bkg
#R divided by G
logratios<-log2(intensities.R[,1:7])-log2(intensities.G[,1:7])
logratios
rownames(logratios)<-derisi$ORF

sgd.orf<-keys(org.Sc.sgd.db,keytype = "ORF")
universe<-intersect(sgd.orf,rownames(logratios))

# logratio of genes that are in universe
universe.logratios<-logratios[which(rownames(logratios) %in% universe),]

jjae1<-universe %in% rownames(logratios)
jjae2<-intersect(rownames(logratios), universe)

a<-apply(universe.logratios,2,mean)

m<-apply(universe.logratios,2,sd)

plot(a,m)
A<-apply(logratios,2,mean)
M<-apply(logratios,2,sd)
plot(a,A)
plot(m,M)
?boxplot

#mean0 and sd0 will be the a[1] and m[1] since time point 0 defines the null distribution

mean0<-a[1]
sd0<-m[1]

calc.zscore<- function(lt1){
  (lt1-mean0)/sd0
}

z1<-sapply(universe.logratios[,1],calc.zscore)
hist(zscores)
hist(zscores7)
qqplot(zscores1,logratios[,1])
# histogram of zscores and logratios at time point 1 seem very similar
z7<-sapply(universe.logratios[,7],calc.zscore)

p1 = pnorm(z1,lower.tail = F)
hist(p1)

p7 = pnorm(z7,lower.tail = F)

disc1<-sum(p7.norm<1e-3) # there are 1104 p-values that are lower than 1e-3
disc1
eval<-length(zscores7)*1e-3#expected values of discoveries that have p val less than 1e-3
eval/disc1
sim.z<-rnorm(length(zscores7))
sim.z.p<-pnorm(-abs(sim.z))
sum(sim.z.p<1e-3) # expected values of discoveries that have p val less than 1e-3 in rnorm

qqplot(zscores7,rnorm(5400,mean=mean(zscores7),sd=sd(zscores7))) # the above result is expected 

disc2<-sum(p7.norm<1e-5) # there are 705 p-values that are lower than 1e-5
eval2<-length(zscores7)*1e-5#expected values of discoveries that have p val less than 1e-5
eval2/disc2
#---------------thomas exercises---------------
x5<-rnorm(5)

zf<-function(n){
  xi<-rnorm(n,0,1)
  (mean(xi)-0)/(1/sqrt(length(xi)))
}

zf(5)
z10k<-replicate(10000,zf(5))
hist(pnorm(z10k))
mean(pnorm(z10k)<0.05)
qqnorm(z10k)
abline(0,1)
  tf<-function(n){
  xi<-rnorm(n,0,1)
  (mean(xi)-0)/(sd(xi)/sqrt(length(xi)))
}

t10k<-replicate(10000,tf(5))
hist(pnorm(t10k))
mean(pnorm(t10k)<0.05)

qqnorm(t10k)

t10k<-replicate(10000,tf(100))
hist(t10k)
qqnorm(t10k)
t10k<-replicate(10000,tf(5))
hist(t10k)
par(mfrow=c(2,2))
