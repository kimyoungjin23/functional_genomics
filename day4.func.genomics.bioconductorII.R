BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
BiocManager::install("airway")
data(airway, package="airway")
se<-airway
assay(se)[56,]
rowRanges(se)[[(rownames(se)[[2]])]] # you dont see rowranges listed when you call for 'se'. but you can see it. why?

dat = assay(se)["ENSG00000120129",]
plot(as.factor(names(dat)),dat)
dextrt = colData(se)$dex

plot(dat, col = dextrt)

library("AnnotationHub")
load("goldenpath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseK562PkRep1.narrowPeak.gz")
ah = AnnotationHub()
mdnase = dnase = query(ah, "goldenpath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseK562PkRep1.narrowPeak.gz")
dnase = query(ah, "goldenpath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseK562PkRep1.narrowPeak.gz")[[1]]
BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
table(seqnames(dnase))
summary(width(dnase))
sum(width(dnase))
sum(width(dnase))
sum(seqlengths(dnase))
sum(width(dnase))/sum(seqlengths(dnase))

# there are three parts to the GRange files, the actual Grange, DataFrame, and Seqinfo

dnase2<-query(ah,"goldenpath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseK562PkRep2.narrowPeak.gz")[[1]]
length(seqnames(dnase2))
summary(width(dnase2))
sum(width(dnase))/sum(seqlengths(dnase))

sum(dnase %over% dnase2)
sum(dnase2 %over% dnase)
sum(ranges(dnase)%over%ranges(dnase2))

w1<-which(dnase %over% dnase2)
w2<-which(ranges(dnase)%over%ranges(dnase2))

dnase[1]
dnase2[2]
dnase2[120043]
ranges(dnase)[1]
ranges(dnase2)[2]
ranges(dnase2)[120043]

grdnase<-GRangesList(
  transcript1 = GRanges(
    seqnames = seqnames(dnase),
    ranges = ranges(dnase)
    ),
  transcript2 = GRanges(
    seqnames = seqnames(dnase2),
    ranges = ranges(dnase2)
  )
)
rgrd<-reduce(unlist(grdnase))

w = width(dnase)
dnase_wide = resize(dnase, width=w+100, fix='center') #make a copy
width(dnase_wide)

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
kg = TxDb.Hsapiens.UCSC.hg19.knownGene
tx= transcripts(kg)
flank(tx, 2000)
proms = promoters(tx)
prom_regions = reduce(proms)
summary(countOverlaps(prom_regions))

findOverlaps(prom_regions, dnase)

prop_proms = sum(width(prom_regions))/sum(seqlengths(prom_regions))
prop_dnase = sum(width(dnase))/sum(seqlengths(prom_regions))
# Iff the dnase and promoter regions are 
# not related, then we would expect this number
# of DNAse overlaps with promoters.
prop_proms * prop_dnase * length(dnase) 
