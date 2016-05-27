require(minfi)
require(minfiData)
load(file="/dcl01/hansen/data/meth_sotos/objects/rgSet.chouf.rda")
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
grSet1.chouf<-grs1
grSet2.chouf<-grs2
grSet3.chouf<-grs3
save(grs1,file="/dcl01/hansen/data/meth_sotos/objects/grSet1.chouf.rda")
save(grs2,file="/dcl01/hansen/data/meth_sotos/objects/grSet2.chouf.rda")
save(grs3,file="/dcl01/hansen/data/meth_sotos/objects/grSet3.chouf.rda")
