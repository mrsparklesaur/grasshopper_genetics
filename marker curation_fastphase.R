# Load libraries ----------------------------------------------------------
library(adegenet)


#File to read
phasefile<-'sanguinipes_fastphase_genotypes.out'
#Read in fastphase, IDs first
fPHASE.IDS <- read.table(phasefile, skip=21, fill=T, sep=" ", comment.char="")[c(TRUE, FALSE, FALSE),]
#Grab IDs and make vector with ID on two rows
IDS <- droplevels(fPHASE.IDS[1:(nrow(fPHASE.IDS)-1),1])
names <- rep(IDS, each = 2)
#Read in genotypes
fPHASE.geno <- droplevels(read.table(phasefile, skip=21, fill=T, sep=" ", comment.char="")[c(FALSE, TRUE, TRUE),])

#Replace brackets (missing data) with -9
fP.mat<-as.matrix(fPHASE.geno)
nadat<-grep("\\[.\\]", fP.mat)
fP.mat[nadat]<- -9
fP.mat<-data.frame(fP.mat)

#Create new data file
dat <- cbind(names,fP.mat)
# write.table(dat,file="sanguinipes.str", row.names=F, col.names=F)

# Read in data in structure format, convert to allele counts --------------
# See https://github.com/ekfchan/evachan.org-Rscripts for scripts
data<-read.structure("sanguinipes.str", n.ind=nlevels(IDS), n.loc=dim(fP.mat)[2], onerowperind=FALSE,
                     col.lab=1, row.marknames=0, col.pop=0, col.others=NULL, NA.char="NA", 
                     ask=FALSE)

geno<-genind2df(data)
geno<-gsub("\"", "", as.matrix(geno)) #these might need to be adjusted; they don't seem to be the same for every file
geno<-gsub("0","", as.matrix(geno)) #same comment as above. Each element should be "NN" where
                                    #N is a nucleotide.
nadat<-grep("-9", geno)
geno[nadat]<-NA

rownames(geno)<-IDS

geno<-t(geno)

source("geno_to_allelecnt.R")

## inputs: matrix of genotypes with rows corresponding to markers and columns to samples
## NA is allowed.
## output: a single matrix of the same size as the input, containing the counts of the
## reference/common allele at each marker.

dfcount<-geno_to_allelecnt(geno)


# Calculate basic SNP statistics ------------------------------------------

# allele frequency (p), MAF (minor allele frequency), 
# MGF (minor genotype frequency), and tests for deviation 
# from HWE (X2 test and Fisher's Exact test)
source("calc_snp_stats.R")

snpstats<-calc_snp_stats(dfcount)

## See also:
# calc_hwe_chisq.R: A script to test for deviation from HWE using Pearson's Chi-Squared test. This test is also incorporated into calc_snp_stats.R.
# calc_wcFst_spop_pairs.R: A script to estimate Fst (theta) values for each pair of sub-populations using the method of Weir & Cockerham 1984 Evolution 38(6): 1358-1370.

## filter by maf
MAFlimit <- which(snpstats$maf>0.02)
mark0<-dfcount[MAFlimit,]

nmark <- ncol(mark0)
missing2 <- apply(mark0,1,function(x){length(which(is.na(x)))/nmark})
mark1b <- mark0[which(missing2 ==0),] #subset data for markers that have <x% missing genotypes
dim(mark1b)

# nsamp<-nrow(mark1b)
# missing1<-apply(mark1b, 2, function(x){length(which(is.na(x)))/nsamp})
# mark2b <- mark1b[,which(missing1 < 0.3)] #subset data for individuals that have <x% missing data
# dim(mark2b)


## filter by HWE
# snpstats2<-calc_snp_stats(mark1b)
# HWEtest<-which(snpstats2$hwe.chisq.p < 0.05)
# mark1<-mark0[HWEtest,]

write.table(mark1, "./sanguinipes_MAF02b.txt", sep="\t")
write.geno(mark1,"sanguinipes_fastphase.geno")

##get locus filter for Jon's code
filter<-read.table("boulderensis_MAF02.txt", sep="\t", as.is=T, check.names=F)
filterloc<-rownames(filter)
filterloc2<-as.numeric(gsub("L","",filterloc))
write.csv(filterloc2,"pellucida_filter.csv")

##set up for arlequin; A,C,T,G data
filter<-read.table("boulderensis_MAF02b.txt", sep="\t", as.is=T, check.names=F)
filterloc<-c("ID","freq",rownames(filter))
colnames(dat)<-c("ID",colnames(geno))
freq<-rep(1, nrow(dat))
dat2<-cbind(dat[1],freq,dat[2:ncol(dat)])

genofilter<-dat2[,colnames(dat2) %in% filterloc]
a<-genofilter[,3:ncol(genofilter)]
b<-apply(a, 1, paste, collapse="")
genofilter2<-cbind(data.frame(genofilter[,1]),
                   data.frame(genofilter[,2]), 
                   data.frame(b))

write.table(genofilter2, "boulderensis_arl2.arp", append=FALSE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

##Set up for arlequin; 0,1,2,3 data
filter<-read.table("sanguinipes_MAF02.txt", sep="\t", as.is=T, check.names=F)
filterloc<-c("ID","freq",rownames(filter))
colnames(dat)<-c("ID",colnames(geno))
freq<-rep(1, nrow(dat))
dat2<-cbind(dat[1],freq,dat[2:ncol(dat)])

genofilter<-dat2[,colnames(dat2) %in% filterloc]
a<-as.matrix(genofilter[,3:ncol(genofilter)])
a2<-gsub("A","0",a)
a2<-gsub("C","1",a2)
a2<-gsub("T","2",a2)
a2<-gsub("G","3",a2)
b<-apply(a2, 1, paste, collapse="")
genofilter2<-cbind(data.frame(genofilter[,1]),
                   data.frame(genofilter[,2]), 
                   data.frame(b))

write.table(genofilter2, "sanguinipes_arl.arp", append=FALSE, quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

##set up for faststructure
dat<-read.table(file="clavatus.str")
dat2<-dat[,2:ncol(dat)] #remove id information. Add this later.
colnames(dat2)<-colnames(geno) #check that geno is in the right orientation
filter<-read.table("clavatus_MAF02.txt", sep="\t", as.is=T, check.names=F)
filterloc<-rownames(filter)
genofilter<-dat2[,colnames(dat2) %in% filterloc]

conv.major.minor<-function(input){
  allele.sums<-lapply(input, table)
  maj.allele<-sapply(allele.sums, function(x) names(x[which.max(x)]))
  min.allele<-sapply(allele.sums, function(x) names(x[which.min(x)]))
  for (i in 1:ncol(input)){
    input[,i]<-gsub(maj.allele[i], "1", input[,i])
    input[,i]<-gsub(min.allele[i], "2", input[,i])
  }
}

output<-conv.major.minor(genofilter)

write.table(output, "clavatus_fstr_12.str", sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)

# Set up for admixture ----------------------------------------------------
filter<-read.table("sanguinipes_MAF02.txt", sep="\t", as.is=T, check.names=F)
filterloc<-rownames(filter)
# make .map file
# 4 columns, L rows (where L is the number of loci)
# Col1=chromosome[integer], Col2=MarkerID[unique string], Col3=Genetic distance[unique float], Col4=physical position[integer]
col1<-rep(1, length(filterloc))
col2<-filterloc
col3<-c(1:length(filterloc))
col4<-c(1:length(filterloc))
map<-cbind(col1, col2, col3, col4)

# make .ped file
# 6 + Lx2 columns; n rows where n is the number of individuals
# Col1=familyID[string], col2=indID[string], col3=fatherID[string], col4=motherID[string],
# col5=sex[int], col6=phenotype[float], col7=locus1 allele1, col8=locus1 allele2.
genofilter<-geno[,colnames(geno) %in% filterloc]
temp<-apply(genofilter, 2, strsplit, "")
pedgen<-t(do.call(mapply, c(cbind, temp)))
colnames(pedgen)<-rep(colnames(genofilter), each=2)
pedgen<-pedgen[match(colnames(filter), rownames(pedgen)),]
metadata<-read.table("./grasshoppers/sanguinipes_sites.txt",sep="\t",comment.char="",as.is=T,check.names=F, header=T)
info<-metadata[match(rownames(pedgen),metadata[,1]),]
pedinfo<-cbind(info$Population, info$Specimen, rep(0, nrow(pedgen)), 
               rep(0, nrow(pedgen)), rep(1, nrow(pedgen)), rep(1, nrow(pedgen)))
ped<-cbind(pedinfo, pedgen)

# write.table(ped, "sanguinipes.ped", col.names=FALSE, row.names=FALSE, quote=FALSE)
# write.table(map, "sanguinipes.map", col.names=FALSE, row.names=FALSE, quote=FALSE)


# # set up for geneland
# filter<-read.table("clavatus_MAF02.txt", sep="\t", as.is=T, check.names=F)
# filterloc<-rownames(filter)
# genofilter<-geno[,colnames(geno) %in% filterloc]
# temp<-apply(genofilter, 2, strsplit, "")
# glnd<-t(do.call(mapply, c(cbind, temp)))
# colnames(glnd)<-rep(colnames(genofilter), each=2)
# 
# write.table(glnd, "clavatus.glnd", col.names=FALSE, row.names=FALSE,
#             quote=FALSE)

