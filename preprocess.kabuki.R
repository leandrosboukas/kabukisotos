load(file="/dcl01/hansen/data/meth_kabuki_bjornsson/objects/rgSet.kabuki.rda")
grSet.kabuki<-preprocessQuantile(rgSet.kabuki)

###keep only the samples from control, KMT2D and KMT2A individuals

grSet.kabuki<-grSet.kabuki[,grSet.kabuki$Phenotype %in% c("normal","KS"),drop=FALSE]
grSet.kabuki<-grSet.kabuki[,grSet.kabuki$KnownGeneMutations %in% c("NoneKnown","MLL2","MLL","MLL4"),drop=FALSE]
save(grSet.kabuki, file="/dcl01/hansen/data/meth_kabuki_bjornsson/objects/grSet.kabuki.rda")