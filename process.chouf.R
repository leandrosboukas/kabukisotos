load(file="/dcl01/hansen/data/meth_sotos/objects/rgSet.chouf.rda")
require(minfi)
require(minfiData)
##now I'm gonna drop the Sotos samples from Hong Kong because they weren't used in the DMP analysis, they were only used later for validation of the methylation signature
hklist<-which(substr(pData(rgSet.chouf)$title,1,2)=="HK") 
hklist<-hklist[-(20:29)]
rgSetinit<-rgSet.chouf[,-hklist,drop=FALSE]
##now keep only control whole blood and sotos whole blood samples
rgSetinit<-rgSetinit[,which(pData(rgSetinit)$group%in%c("Control","Sotos")),drop=FALSE]
rgSetinit<-rgSetinit[,which(pData(rgSetinit)$tissue%in%c("whole blood")),drop=FALSE]
##now create three different rgSets because they did 3 comparisons of the 53 controls with 17 sotos patients, each time containing 16 sotos and one of the three family members who all had sotos
rg1<-rgSetinit[,-which(substr(rgSetinit$title,1,7)%in%c("DL50450","DL50452")),drop=FALSE]
rg2<-rgSetinit[,-which(substr(rgSetinit$title,1,7)%in%c("DL50448","DL50452")),drop=FALSE]
rg3<-rgSetinit[,-which(substr(rgSetinit$title,1,7)%in%c("DL50448","DL50450")),drop=FALSE]
##now process the 3 rgsets using the illumina normalization because that's what they did, and then find the DMPs
M1<-preprocessIllumina(rg1,bg.correct=TRUE,normalize="controls")
M2<-preprocessIllumina(rg2,bg.correct=TRUE,normalize="controls")
M3<-preprocessIllumina(rg3,bg.correct=TRUE,normalize="controls")
rset1<-ratioConvert(M1, what = "both", keepCN = TRUE)
grs1<- mapToGenome(rset1)
rset2<-ratioConvert(M2, what = "both", keepCN = TRUE)
grs2<-mapToGenome(rset2)
rset3<-ratioConvert(M3, what = "both", keepCN = TRUE)
grs3<-mapToGenome(rset3)
beta1<-getBeta(grs1)
beta2<-getBeta(grs2)
beta3<-getBeta(grs3)
group1<-pData(grs1)$group
group2<-pData(grs2)$group
group3<-pData(grs3)$group
##now do the 3 different comparisons to get the unified list
dmp1<-dmpFinder(beta1,group1,type="categorical",qCutoff=0.05)
dmp2<-dmpFinder(beta2,group2,type="categorical",qCutoff=0.05)
dmp3<-dmpFinder(beta3,group3,type="categorical",qCutoff=0.05)
dmplistsot<-Reduce(intersect,list(rownames(dmp1),rownames(dmp2),rownames(dmp3)))
