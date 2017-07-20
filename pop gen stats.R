library("adegenet")
library("hierfstat")
library("pegas")
library("ade4")

source("C:/Users/Rachel/Documents/Work files/R/GrasshopperGenetics/basic.stats2.R")
source("C:/Users/Rachel/Documents/Work files/R/GrasshopperGenetics/pop gen helper functions.R")
source("C:/Users/Rachel/Documents/Work files/R/GrasshopperGenetics/pairwise.neifst2.R")

# load data and convert to genind format ----------------------------------

data<-read.table("./sanguinipes_MAF02.txt", sep="\t", as.is=T, check.names=F, header=T)
#transpose data
data<-t(data)

metadata <- read.table("./grasshoppers/sanguinipes_sites.txt",sep="\t",comment.char="",as.is=T,check.names=F, header=T)
info <- metadata[match(rownames(data),metadata[,1]),]

# datatemp<-data[which(rownames(data) %in% info[,1]),]
# info<-info[match(rownames(datatemp), info[,1]),]

ind<-as.character(info$Specimen)
population<-as.character(info$Population)

# colnames(data)<-gsub("\\.","_", colnames(data))

#convert to genind
# gendat<-df2genind(data, ploidy=2, ind.names=ind, pop=population, sep="",NA.char=9)

data3split<-split(data.frame(data), population)

#convert to "loci"
data2<-data.frame(data) #adegenet
# loci<-as.loci(data2) #pegas

#write STRUCTURE file
# write.loci(data2, file="sanguinipes_data",loci.sep="\t",quote=FALSE,
# allele.sep="\t", na="-9", col.names=FALSE)

data3<-cbind(data.frame(population), data2)

#fill missing values with mode
# impute.mode <- function(x) {ux <- unique(x);ux[which.max(tabulate(match(x, ux)))]}
# dataM <- apply(data,1,function(x){ix <- which(is.na(x));x[ix] <- impute.mode(x); return(x)})
# dataM<-t(data.frame(dataM))
# colnames(dataM)<-gsub("\\.","_", colnames(dataM))
# gendatM<-df2genind(dataM, ploidy=2, ind.names=ind, pop=population, sep="")

# Basic statistics --------------------------------------------------------

# ## genetic diversity
# div<-summary(gendat)
# div2<-summary(loci)
# #population-specific estimates of Hexp
# genpop<-lapply(data3split, df2genind, ploidy=2, sep="",NA.char=9)
# gendiv<-lapply(genpop, summary)
# genhe<-sapply(gendiv, function(x) mean(x$Hexp))
# 
# #heterozygosity
# het<-sapply(div2, function(x) H(x$allele))
# mean(het) #this is nearly (but not quite) the same as mean(div$Hexp)
# 
# #mean heterozygosity
# mean(div$Hexp)

# #calculate genetic diversity
# thetaH <- sapply(div2, function(x) theta.h(x$allele))
# mean(thetaH)
# thetaK <- sapply(div2, function(x) theta.k(x$allele))
# mean(thetaK)

bstat<-basic.stats2(data3)
bstat$overall

popHo<-apply(bstat$Ho, 2, mean, na.rm=TRUE)
print(popHo, digits=3)
popHs<-apply(bstat$Hs, 2, mean, na.rm=TRUE)
print(popHs, digits=3)
popFis<-apply(bstat$Fis, 2, mean, na.rm=TRUE)
print(popFis, digits=3)


# My own way of calculating observed heterozygosity -----------------------
mean(calc.Hobs(data2))
print(sapply(calc.Hobs.pop(data2, population), mean), digits=4)

# My own way of calculating expected heterozygosity -----------------------
mean(calc.Hexp(data2))
print(sapply(calc.Hexp.pop(data2, population), mean), digits=4)

# Count number of polymorphic loci in each population
sapply(data3split, count.polymorphic)

## test for HWE
# hwtest<-hw.test(gendat, B = 0) # doesn't work with raw data - just producing NA

# Population differentiation ----------------------------------------------
# fstat(gendat)
# gsum<-gstat.randtest(gendat, nsim=99) #doesn't really work with missing data

# pairwise Fst
# genet.dist(gendat, method="Nei87")
pairwise.neifst2(data3) #Fst according to Nei 1987
# pairwise.WCfst(data3) #Fst according to Weir & Cockerham 1984

# # DAPC --------------------------------------------------------------------
# library(RColorBrewer)
# bg<-brewer.pal(5, "Dark2")
# 
# mb.dapc<-dapc(gendat, n.pca=30, pca.select="nbEig", n.da=3)
# summary(mb.dapc) #look esp at assign.prop, assign.per.pop, and post.grp.size
# 
# scatter(mb.dapc, col=bg, scree.da=FALSE, cstar=FALSE)
# 
# ac.dapc<-dapc(gendat, n.pca=35, pca.select="nbEig", n.da=3)
# summary(ac.dapc)
# scatter(ac.dapc, col=bg, scree.da=FALSE, cstar=FALSE, label=NULL)
# 
# ms.dapc<-dapc(gendat, n.pca=35, pca.select="nbEig", n.da=4)
# summary(ms.dapc)
# scatter(ms.dapc, col=bg, scree.da=FALSE, cstar=FALSE, label=NULL)
# 
# cp.dapc<-dapc(gendat, n.pca=30, pca.select="nbEig", n.da=2)
# summary(cp.dapc)
# scatter(cp.dapc, col=bg, scree.da=FALSE, cstar=FALSE)
# 
# 
# # cross-validation
# test<-xvalDapc(data, grp=population, n.pca.max=80)
