dds.n <- makeExampleDESeqDataSet(n=1000, m=4)
dds.n <- estimateSizeFactors(dds.n)
sizeFactors(dds.n)

dds.n <- estimateSizeFactors(dds.n, controlGenes=1:200)

m <- matrix(runif(1000 * 4, .5, 1.5), ncol=4)
dds.n <- estimateSizeFactors(dds.n, normMatrix=m)
normalizationFactors(dds.n)[1:3,]

geoMeans <- exp(rowMeans(log(counts(dds.n))))
dds.n <- estimateSizeFactors(dds.n,geoMeans=geoMeans)
sizeFactors(dds.n)