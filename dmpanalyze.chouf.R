load(file="/dcl01/hansen/data/meth_sotos/objects/grSet1.chouf.rda")
load(file="/dcl01/hansen/data/meth_sotos/objects/grSet2.chouf.rda")
load(file="/dcl01/hansen/data/meth_sotos/objects/grSet3.chouf.rda")
beta1<-getBeta(grSet1.chouf)
beta2<-getBeta(grSet2.chouf)
beta3<-getBeta(grSet3.chouf)
group1<-pData(grSet1.chouf)$group
group2<-pData(grSet2.chouf)$group
group3<-pData(grSet3.chouf)$group
##now do the 3 different comparisons and then get the unified list
dmp1<-dmpFinder(beta1,group1,type="categorical",qCutoff=0.05)
dmp2<-dmpFinder(beta2,group2,type="categorical",qCutoff=0.05)
dmp3<-dmpFinder(beta3,group3,type="categorical",qCutoff=0.05)
dmplistsot<-Reduce(intersect,list(rownames(dmp1),rownames(dmp2),rownames(dmp3)))

##now do the same dmp analysis using the mann whitney test with bonferroni correction because that's what they did
mwcont<-which(grSet1.chouf$group=="Control")
mwsot<-which(grSet1.chouf$group=="Sotos")
##define functions that perform the mann whitney test for each row of the beta matrix and return the p value
mwfunc1<-function(numeric){return(wilcox.test(beta1[numeric,mwcont],beta1[numeric,mwsot])[3])}
mwfunc2<-function(numeric){return(wilcox.test(beta2[numeric,mwcont],beta2[numeric,mwsot])[3])}
mwfunc3<-function(numeric){return(wilcox.test(beta3[numeric,mwcont],beta3[numeric,mwsot])[3])}
dmpmw1<-sapply(1:485512,mwfunc1)
dmpmw2<-sapply(1:485512,mwfunc2)
dmpmw3<-sapply(1:485512,mwfunc3)
##do the bonferroni correction
dmpmwbonf1<-p.adjust(dmpmw1,method="bonferroni")
dmpmwbonf2<-p.adjust(dmpmw2,method="bonferroni")
dmpmwbonf3<-p.adjust(dmpmw3,method="bonferroni")
dmpmw1<-which(dmpmwbonf1<0.05)
dmpmw2<-which(dmpmwbonf2<0.05)
dmpmw3<-which(dmpmwbonf3<0.05)
##get cpg identities
givecpg1<-function(numeric){return(rownames(beta1)[dmpmw1[numeric]])}
givecpg2<-function(numeric){return(rownames(beta2)[dmpmw2[numeric]])}
givecpg3<-function(numeric){return(rownames(beta3)[dmpmw3[numeric]])}
dmpsot1<-sapply(1:length(dmpmw1),givecpg1)
dmpsot2<-sapply(1:length(dmpmw2),givecpg2)
dmpsot3<-sapply(1:length(dmpmw3),givecpg3)
mwdmplistsot<-Reduce(intersect,list(dmpsot1,dmpsot2,dmpsot3))
##verify the dmps identified by mann whitney are contained in the previous list I got after dmpFinder
length(which(mwdmplistsot%in%dmplistsot))