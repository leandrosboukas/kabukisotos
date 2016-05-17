load(file="/dcl01/hansen/data/meth_sotos/objects/rgSet.sotos.rda")
load(file="/dcl01/hansen/data/meth_sotos/objects/grSet.sotos.rda")
gr<-granges(grSet.sotos)
an<-getAnnotation(grSet.sotos)
beta <- getBeta(grSet.sotos)
disgroup  <- pData(grSet.sotos)$group
dmp <- dmpFinder(beta, pheno = disgroup, type = "categorical", qCutoff=0.05)
head(dmp)
tail(dmp)
##import Hans' Kabuki data
library(gdata)
dmpkab<-read.xls("/home/other/lboukas/DMPs Kabuki type 1 vs controls Hans' data.xls")
head(dmpkab)
##function for cpg identity from chr position on chromosome
cpgident<-function(numeric){for (i in 1:485512){if (an$pos[i]==dmpkab[numeric,2]){return(rownames(an)[i])}}}
##create a list of the cpg identities of the Kabuki vs controls dmp's
cpglist<-vector()
for (j in 1:57){cpglist<-append(cpglist,cpgident(j))}
cpglist
##functions that output chromosome number and cpg position on chromosome to verify that the cpglist from above is correct
cpgchr<-function(character){for (i in 1:485512){if (rownames(an)[i]==character)return(an$chr[i])}}
cpgpos<-function(character){for (i in 1:485512){if (rownames(an)[i]==character)return(an$pos[i])}}
cpgchr(cpglist[1])
cpgpos(cpglist[1])
##function that outputs row of sotos dmp table from cpg identity of Kabuki dmp (or NULL if it is not a shared dmp)
cpgqsot<-function(character){for (i in 1:310089){if (rownames(dmp)[i]==character)return(dmp[i,])}}
##create table for shared dmps
dmpshared<-data.frame()
for (j in 1:57){dmpshared<-rbind(dmpshared,cpgqsot(cpglist[j]))}
dim(dmpshared)
head(dmpshared)
##function that outputs the difference of average beta values (sotos - kabuki) for a given cpg
cpgbetadif<-function(character){for (i in 1:485512){if (rownames(an)[i]==character)return((sum(beta[i,54:91])/38)-(sum(beta[i,1:53])/53))}}
##create list of beta differences and add it as a column to dmpshared
deltabetalist<-vector()
for (j in 1:52){deltabetalist<-append(deltabetalist,cpgbetadif(rownames(dmpshared)[j]))}
deltabetalist
dmpshared<-cbind(dmpshared,deltabetalist)
head(dmpshared)
######the rest is just manipulations so that I can create a combined table with the Kabuki and Sotos dmps
cpgposlist<-vector()
for (j in 1:52){cpgposlist<-append(cpgposlist,cpgpos(rownames(dmpshared)[j]))}
head(cpgposlist)
#remove dmps that are not shared dmps with Sotos from dmpkab
for (j in 1:57){if ((dmpkab$Position[j] %in% cpgposlist)!=TRUE){print(j)}}
dmpkab2<-dmpkab[-c(19,20,30,31,42),]
head(dmpkab2)
#verify, the following should not output anything
for (j in 1:52){if (cpgpos(rownames(dmpshared)[j])!=dmpkab2$Position[j]){print(j)}}
#add the cpg identity column
dmpkab2<-cbind(dmpkab2,rownames(dmpshared))
#reorder columns
dmpkab2<-dmpkab2[,c(1,2,8,3,4,5,6,7)]
colnames(dmpkab2)[3]<-"CpG identity"
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
