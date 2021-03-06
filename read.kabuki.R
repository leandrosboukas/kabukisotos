library(minfi)
targets <- read.metharray.sheet("/dcl01/hansen/data/meth_kabuki_bjornsson", "Annotation_v2.csv")
targets
rgSet <- read.metharray.exp(targets = targets, verbose = TRUE)

rgSet.kabuki <- rgSet
save(rgSet.kabuki, file="/dcl01/hansen/data/meth_kabuki_bjornsson/objects/rgSet.kabuki.rda")