# course material website: https://seandavi.github.io/ITR/

rnorm(100)
my_rn<-rnorm(100)
mean(my_rn)
sd(my_rn)
my_rn+5
mean(my_rn+5)
speakers<-sample(c(1:25))
speakers
speakers<-sample(c(1:25),3)
speakers
speakers<-sample(c(1:25),3)
speakers
speakers<-sample(c(1:25),3)
speakers<-sample(c(1:25))
day1<-speakers[1:13]
day2<-speakers[14:25]
day1
rep("A",10)
my_rn[sd(my_rn)>1]
my_rn[c(sd(my_rn)>1)]
sd(my_rn)>1
my_rn>(mean(my_rn)-sd(my_rn))
my_rn[my_rn>(mean(my_rn)-sd(my_rn))]
my_rn.onesd<-my_rn[my_rn>(mean(my_rn)-sd(my_rn))]
hist(my_rn)
abline(v= c((mean(my_rn)-sd(my_rn)),(mean(my_rn)+sd(my_rn)),mean(my_rn)),col= c(2,2,3))

install.packages("tidyverse")
library(tidyverse)

ydat <- read.csv('https://raw.githubusercontent.com/bioconnector/workshops/master/data/brauer2007_tidy.csv')
leucine.bp<-filter(ydat,nutrient=="Leucine",bp=="leucine biosynthesis")
head(leucine.bp)
dim(leucine.bp)
quantile(ydat$expression, probs=.99)
quantile(ydat$expression, 2.07)
hist(ydat$expression)
exp.99<-which(ydat$expression>2.07)
view(ydat$bp[exp.99]) # functions of genes expressed at 99th percentile expression


nogo <- select(ydat, -bp, -mf)
nogo

# we could filter this new dataset
filter(nogo, symbol=="LEU1" & rate==.05)

# Notice how the original data is unchanged - still have all 7 columns
ydat

mutate (nogo, signal = 2^expression)
new.nogo<-mutate (nogo, signal = 2^expression, sqrtsg = sqrt(signal))
arrange(new.nogo, signal)

summarize(group_by(ydat, symbol), meanexp= mean(expression))

filter(ydat, symbol == "ADH2"& rate==0.05) %>% select(expression)

filter(ydat, rate == 0.05, bp == "response to stress") %>% group_by(nutrient) %>% 
  summarize(meanexp= mean(expression))

summarize(ydat, n_distinct(mf))
group_by(ydat, bp) %>% summarize(frequency=n_distinct(symbol)) %>% arrange (desc(frequency))
filter(ydat, bp!="biological process unknown" & mf=="molecular function unknown")%>% 
  select(symbol, bp, mf) %>%
  distinct()

ydat %>% 
  filter(bp!="biological process unknown" & mf=="molecular function unknown") %>% 
  select(symbol, bp, mf) %>% 
  summarize(bp=n_distinct(bp), genes = n_distinct(symbol), molecular_function = n_distinct(mf))

filter(ydat, rate == 0.05) %>% group_by(bp) %>% summarize(meanexp=mean(expression)) %>%
  mutate(meanexp.2=round(meanexp,2))
  arrange(desc(meanexp))

r100<-rnorm(100)
mean(r100)
var(r100)
sd(r100)
abs(sd(r100)^2-var(r100))
r10k<-rnorm(10000)
summarize(r10k)
summary(r10k)
plot.ecdf(r10k)
plot(pnorm(seq(-5,5,along.with = r10k)))
r10k.s2<-rnorm(1e4,0,2)
par(mfrow=c(1,1))
plot.ecdf(r100)
plot.ecdf(r10k, add=TRUE)

abline(v=1)
abline(h=0.8, cex=1)
plot.ecdf(r10k.s2, col = "blue", add=TRUE)
r10k.m1<-rnorm(1e4,1,1)
plot.ecdf(r10k.m1,  col = "red", add= TRUE, verticals = TRUE,col.vert = "green")

r10K.lt0<-r10k<0
head(r10K.lt0)
table(r10K.lt0)
pnorm(0)
r10K.gt1<-r10k>1
r10K.gt1
sum(r10K.gt1)
table(r10K.gt1)
1-pnorm(1)
sum(r10k.s2>10)
1-pnorm(10,0,2)

abline(h=.5)

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
logratios[ACR1,]
hours<-seq(0,12, by=2)
plot(hours,logratios[ACR1,],type = "l")
#downregulated gene
SAM1<-which(derisi$Name=="SAM1")
plot(hours,logratios[SAM1,], type = "l",add=TRUE)

order(logratios[,7]-logratios[,1], decreasing = TRUE)
(logratios[,7]-logratios[,1])[5]
time.logratios<-logratios[,7]-logratios[,1]
derisi$Name[which(time.logratios==max(time.logratios))] #most upregulated gene
derisi$Name[which(time.logratios==min(time.logratios))] #most downregulated gene
derisi$Name[2340:2359]
time.means<-apply(logratios, 2, mean)
plot(time.means)

hist(logratios[,1])

#MA plot x = A: mean log expression. y= mean log difference
M=logratios[,1]
A=(log2(intensities.R[,1])+log2(intensities.G[,1]))/2
plot(A,M)
cor(A,M)
cor(logratios[,1],logratios[,2])
plot(logratios[,1],logratios[,2])

plot.ecdf(logratios[,1])
plot.ecdf(logratios[,7],add=TRUE,col="blue")

cor(logratios[,1:7])
round(var(logratios),2)
apply(logratios,2,var)
plot(hours,apply(logratios,2,var))

