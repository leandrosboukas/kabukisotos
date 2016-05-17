library(minfi)
targets <- read.metharray.sheet("/dcl01/hansen/data/meth_kabuki_bjornsson", "Annotation_v2.csv")
targets
rgSet <- read.metharray.exp(targets = targets, verbose = TRUE)

## CHECK AND SAVE
