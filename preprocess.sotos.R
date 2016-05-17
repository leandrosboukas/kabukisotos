load(file="/dcl01/hansen/data/meth_sotos/objects/rgSet.sotos.rda")
grSet.sotos<-preprocessQuantile(rgSet.sotos)

###keep only whole blood samples from Control and Sotos individuals

grSet.sotos<-grSet.sotos[,grSet.sotos$group %in% c("Control","Sotos"),drop=FALSE]
grSet.sotos<-grSet.sotos[,grSet.sotos$tissue %in% c("whole blood"),drop=FALSE]
save(grSet.sotos, file="/dcl01/hansen/data/meth_sotos/objects/grSet.sotos.rda")