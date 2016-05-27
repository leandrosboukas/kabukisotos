require(minfi)
require(minfiData)
library(GEOquery)
idatFiles <- list.files("GSE74432/idat", pattern = "idat.gz$", full = TRUE)
rgSet <- read.450k.exp("GSE74432/idat",verbose=TRUE)
rgSet
pData(rgSet)
head(sampleNames(rgSet))
geoMat <- getGEO('GSE74432',GSEMatrix=TRUE)
pD <- pData(geoMat[[1]])
sampleNames(rgSet) <- sub(".*_7", "7", sampleNames(rgSet))
sampleNames(rgSet) <- sub(".*_8", "8", sampleNames(rgSet))
sampleNames(rgSet) <- sub(".*_9", "9", sampleNames(rgSet))
sampleNames(rgSet) <- sub(".*_3", "3", sampleNames(rgSet))
names(pD)[c(11,12)] <- c("group", "tissue")
pD$group <- sub("^disease state: ", "", pD$group)
pD$tissue <- sub("^tissue: ", "", pD$tissue)
pD<-cbind(pD,sampleNames(rgSet))
pData(rgSet)<-pD
colnames(pData(rgSet))[34]<-"sampleNames"
rgSet.chouf<-rgSet

## New
pData(rgSet.chouf)$Array <- sub(".*_", "", pData(rgSet.chouf)$sampleNames)
pData(rgSet.chouf)$Slide <- sub("_.*", "", pData(rgSet.chouf)$sampleNames)
## Remove all the stupid columns we will never use
## Rename sampleNames
## Clean up gender



save(rgSet.chouf, file="/dcl01/hansen/data/meth_sotos/objects/rgSet.chouf.rda")
