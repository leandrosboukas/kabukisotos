library(minfi)
library(GEOquery)
idatFiles <- list.files("/dcl01/hansen/data/meth_sotos/GSE74432/idat", pattern = "idat.gz$", full = TRUE)
rgSet <- read.450k.exp("/dcl01/hansen/data/meth_sotos/GSE74432/idat", verbose=TRUE)
rgSet
pData(rgSet)
head(sampleNames(rgSet))
geoMat <- getGEO('GSE74432', GSEMatrix=TRUE)
pD.all <- pData(geoMat[[1]])
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
head(pD)
names(pD)[c(3,4)] <- c("group", "tissue")
pD$group <- sub("^disease state: ", "", pD$group)
pD$tissue <- sub("^tissue: ", "", pD$tissue)
sampleNames(rgSet) <- sub(".*_7", "7", sampleNames(rgSet))
sampleNames(rgSet) <- sub(".*_8", "8", sampleNames(rgSet))
sampleNames(rgSet) <- sub(".*_9", "9", sampleNames(rgSet))
sampleNames(rgSet) <- sub(".*_3", "3", sampleNames(rgSet))
head(sampleNames(rgSet))
pD$title <- sampleNames(rgSet)
rownames(pD) <- pD$title
pD <- pD[sampleNames(rgSet),]
pData(rgSet) <- pD  
rgSet

rgSet.sotos <- rgSet
save(rgSet.sotos, file="/dcl01/hansen/data/meth_sotos/objects/rgSet.sotos.rda")



