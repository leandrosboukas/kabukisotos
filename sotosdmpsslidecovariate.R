require(minfi)
require(minfiData)
library(siggenes)

load(file="/dcl01/hansen/data/meth_sotos/objects/rgSet.sotos.rda")
pData(rgSet.sotos)$Array <- sub(".*_", "", sampleNames(rgSet.sotos))
pData(rgSet.sotos)$Slide <- sub("_.*", "", sampleNames(rgSet.sotos))
rgSet.to.use<-rgSet.sotos[,rgSet.sotos$group%in%c("Control","Sotos"),drop=FALSE]
rgSet.to.use<-rgSet.to.use[,rgSet.to.use$tissue%in%c("whole blood"),drop=FALSE]

grset<-preprocessQuantile(rgSet.to.use)
gr<-granges(grset)
an<-getAnnotation(grset)
betas<-getBeta(grset)
ph1<-grset$group
ph2<-grset$Slide
#save(grset,file="home/other/lboukas/sotos.grset.slide.covariate.rda")

###define function to perform rergression with group and slide as predictor variables and beta value as response variable
regrcpg<-function(numeric){return(lm(betas[numeric,]~as.factor(ph1)+as.factor(ph2),data=pData(grset)))}
#define functions that output the intercept, and the p value for the coefficient of the group variable
intfun<-function(numeric){return(summary(regrcpg(numeric))$coefficients[1,1])}
pvalfun<-function(numeric){return(summary(regrcpg(numeric))$coefficients[2,4])}
#get the list of intercepts and p values (the original dmpFinder function also outputs the F statistic, can get that with summary(regrcpg)[10] but I don't think I need it)
interceptlist<-sapply(1:dim(betas)[1],intfun)
pvallist<-sapply(1:dim(betas)[1],pvalfun)
#get the list of q values (as done in the siggenes package)
p0 <- pi0.est(pvallist)$p0
qvallist<-qvalue.cal(pvallist,p0)
#get the table of DMPs (if I did regression with only the group as predictor variable, this table would be the same as the output from dmpFinder (minus the F statistic column))
dmpstab <- data.frame(intercept=interceptlist,pval=pvallist,qval=qvallist)
rownames(dmpstab)<-rownames(betas)
sotosdmps<-dmpstab[order(dmpstab$qval),]
sotosdmps<-sotosdmps[1:length(which(dmpstab$qval<0.05)),]
##function that outputs the difference of average beta values (sotos - controls) for a given cpg
cpgbetadif<-function(numeric){return((sum(betas[numeric,54:91])/38)-(sum(betas[numeric,1:53])/53))}
#create list of deltabeta differences for the sotos dmps only
dbhelplist2<-which(rownames(betas)%in%rownames(sotosdmps))
sotdeltabetalist<-sapply(dbhelplist2,cpgbetadif)
###save some objects
#save(interceptlist,file="/home/other/lboukas/intlist.rda")
#save(pvallist,file="/home/other/lboukas/pvallist.rda")
#save(sotosdmps,file="/home/other/lboukas/sotosdmpstable.rda")
#save(deltabetalist,file="/home/other/lboukas/deltabetalistsot.rda")
##import Hans' Kabuki data
library(gdata)
dmpkab<-read.xls("/home/other/lboukas/DMPs Kabuki type 1 vs controls Hans' data.xls")
##function for cpg identity from chr position on chromosome
cpgident<-function(numeric){for (i in 1:485512){if (an$pos[i]==dmpkab[numeric,2]){return(rownames(an)[i])}}}
##create a list of the cpg identities of the Kabuki vs controls dmp's
cpglist<-sapply(1:length(rownames(dmpkab)),cpgident)
cpglist
##functions that output chromosome number and cpg position on chromosome to verify that the cpglist from above is correct
cpgchr<-function(character){for (i in 1:485512){if (rownames(an)[i]==character)return(an$chr[i])}}
cpgpos<-function(character){for (i in 1:485512){if (rownames(an)[i]==character)return(an$pos[i])}}
cpgchr(cpglist[1])
cpgpos(cpglist[1])
#change the rownames of dmpkab list so that each rowname is the cpg identity
rownames(dmpkab)<-cpglist
##function that outputs row of sotos dmp table from cpg identity of Kabuki dmp (or NULL if it is not a shared dmp)
cpgqsot<-function(character){for (i in 1:length(rownames(sotosdmps))){if (rownames(sotosdmps)[i]==character)return(sotosdmps[i,])}}
##create table for shared dmps
dmpshared<-data.frame()
for (j in 1:57){dmpshared<-rbind(dmpshared,cpgqsot(cpglist[j]))}
dim(dmpshared)
#create list of deltabeta differences for the shared dmps and add it as a column to the shared dmp table
dbhelplist<-which(rownames(betas)%in%rownames(dmpshared))
shdeltabetalist<-sapply(dbhelplist,cpgbetadif)
dmpshared<-cbind(dmpshared,deltabetalist)


######the rest is just manipulations so that I can create a combined table with the Kabuki and Sotos dmps
cpgposlist<-vector()
for (j in 1:dim(dmpshared)[1]){cpgposlist<-append(cpgposlist,cpgpos(rownames(dmpshared)[j]))}
#remove dmps that are not shared dmps with Sotos from dmpkab
dmpkab2<-dmpkab[which(rownames(dmpkab)%in%rownames(sotosdmps)),]
head(dmpkab2)
#verify, the following should not output anything
for (j in 1:52){if (cpgpos(rownames(dmpshared)[j])!=dmpkab2$Position[j]){print(j)}}
colnames(dmpkab2)[5]<-"deltaBeta Kabuki vs Controls"
colnames(dmpkab2)[6]<-"p val Kab"
colnames(dmpkab2)[7]<-"q val Kab"
##add the deltabeta Sotos column
dmpkab2<-cbind(dmpkab2,dmpshared$deltabetalist)
#add the Sotos p and q value columns
dmpkab2<-cbind(dmpkab2,dmpshared$pval)
dmpkab2<-cbind(dmpkab2,dmpshared$qval)
colnames(dmpkab2)[9]<-"deltaBeta Sotos vs Controls"
colnames(dmpkab2)[10]<-"p val Sot"
colnames(dmpkab2)[11]<-"q val Sot"
head(dmpkab2)
library(xlsx)
write.xlsx(dmpkab2, "/home/other/lboukas/shareddmpskabukisotos.xlsx")







