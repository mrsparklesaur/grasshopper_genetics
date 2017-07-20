###########################################################################
# ABC analysis ------------------------------------------------------------
###########################################################################

# Functions
extract.files<-function(tarfiles){
    dir.create("params_extract")
    dir.create("sumstats_extract")
    for (i in 1:length(tarfiles)){
        untar(tarfiles[[i]], exdir="temp_extract")
        # check size of param files
        checkfiles<-list.files(path="./temp_extract", pattern=".*params", full.names=TRUE)
        size<-file.info(checkfiles)
        keepfiles<-gsub(".*temp_extract/(.*).params", "\\1", rownames(size)[which(size$size!=0)])
        statfiles<-list.files(path="./temp_extract", pattern=".sumstats", full.names=TRUE)
        statfiles2<-gsub(".*temp_extract/(.*).sumstats", "\\1", statfiles)
        loading<-statfiles2[which(statfiles2 %in% keepfiles)]
        file.copy(paste("./temp_extract/", loading, ".params", sep=""), "./params_extract", copy.mode=TRUE)
        file.copy(paste("./temp_extract/", loading, ".sumstats", sep=""), "./sumstats_extract", copy.mode=TRUE)
        unlink("temp_extract", recursive=TRUE)
    }
}

# Packages
library(data.table)
library(abc)

# Read in files -----------------------------------------------------------

# Unpack tar.gz files
tarfiles<-list.files(pattern="tar.gz")
extract.files(tarfiles)

# Read in .params files
paramfiles<-list.files(path="./params_extract",pattern="*.params", full.names=TRUE)

paramlist<-list()
for (i in 1:length(paramfiles)){
  paramlist[[i]]<-read.table(paramfiles[i], header=TRUE, fill=TRUE)[,1:5]
}

paramnames<-gsub(".*\\d+_([a-z]*)_\\d+.params","\\1", paramfiles)
names(paramlist)<-substr(paramnames, 3, 4)

# Read in .sumstats files
statfiles<-list.files(path="./sumstats_extract",pattern="*.sumstats", full.names=TRUE)
# Check that files have a file size
size<-file.info(statfiles)
size<-size$size
statfiles<-statfiles[which(size!=0)]
paramlist<-paramlist[which(size!=0)]

statlist<-list()
for (i in 1:length(statfiles)){
  statlist[[i]]<-read.table(statfiles[i], header=TRUE, fill=TRUE)
}

statnames<-gsub(".*\\d+_([a-z]*)_\\d+.sumstats","\\1", statfiles)
names(statlist)<-substr(statnames, 3, 4)

# Check that the length of the params and sumstats lists match
check<-vector()
for (i in 1:length(statlist)){
  a<-nrow(statlist[[i]])
  b<-nrow(paramlist[[i]])
  if (a == b){
    check[[i]]<-1
  } else {
    check[[i]]<-0
  }
}

paramlist<-paramlist[which(check==1)]
statlist<-statlist[which(check==1)]

# Create params and sumstats data frames
par.sim<-rbindlist(paramlist)
namelength<-sapply(paramlist, nrow)
models<-unname(unlist(mapply(rep, names(paramlist), each=namelength)))
par.sim$models<-models
# write.table(par.sim, "cp.paramsim", row.names=F, quote=F, sep="\t")

stat.sim<-rbindlist(statlist)
stat.sim$models<-models
# Remove summary stats that have zero variance
novar<-which(sapply(stat.sim[,1:(ncol(stat.sim)-1)], var)==0)
stat.sim<-as.data.frame(stat.sim)[,-novar]
# write.table(stat.sim2, "cp.statsim", row.names=F, quote=F, sep="\t")

# ABC ---------------------------------------------------------------------

# par.sim<-read.table("cp.paramsim", header=TRUE, sep="\t")
# stat.sim<-read.table("cp.statsim", header=TRUE, sep="\t")

# Read in observed summary statistics (as a vector)
stats.obs<-read.table("cpoutput.sumstats", header=TRUE, fill=TRUE)
stats.obs<-stats.obs[,-novar]

# check number of models and truncate to 2 million
statsplit<-split(stat.sim, f=stat.sim$models)
modlength<-sapply(statsplit, function(x) nrow(x))
stattrunc<-lapply(statsplit, function(x) x[1:2000000,])
stat.sim2<-rbindlist(stattrunc)
t<-ncol(stat.sim2)-1

# truncate param files
paramsplit<-split(par.sim, f=par.sim$models)
paramtrunc<-lapply(paramsplit, function(x) x[1:2000000,])
par.sim2<-rbindlist(paramtrunc)

# Vector of model indices
models<-as.vector(stat.sim2$models)
mods<-unique(models)

# Look at summary stats
png("summary_stats.png", width=600, height=600)
par(mfcol = c(2,3), mar=c(5,3,4,0.5))
boxplot(stat.sim2$mean_Pi~models, main="Mean nucleotide diversity")
boxplot(stat.sim2$sd_Pi~models, main="SD nucleotide diversity")
boxplot(stat.sim2$mean_H~models, main="Mean heterozygosity")
boxplot(stat.sim2$sd_H~models, main="SD heterozygosity")
boxplot(stat.sim2$tot_H~models, main="Total Heterozygosity")
boxplot(stat.sim2$mean_S~models, main="Mean S")

# Cross-validation: can ABC distinguish between the models (using multinomial logistic regression)?
# change nval to 1000 for real run
# I think these tolerance values are in %, so 0.1 keeps 0.1% of the accepted simulations
# Tolerance values keep (x) of the accepted simulations, so 0.01 will keep 1% of the accepted simulations
cv.modsel<-cv4postpr(models, stat.sim2[,1:t], nval=100, tols=c(0.0005, 0.001, 0.002, 0.005, 0.01, 0.02), method="neuralnet")
cv.modsel2<-cv4postpr(models, stat.sim2[,1:t], nval=100, tols=c(0.0005, 0.001, 0.002, 0.005, 0.01, 0.02), method="mnlogistic")
s<-summary(cv.modsel)
plot(cv.modsel, names.arg=mods)

# # Cross-validation: can ABC distinguish between the models (using rejection method)?
# cv.modsel.rej<-cv4postpr(models, stat.sim, nval=10, tols=c(0.01, 0.05, 0.1, 0.2), method="rejection")
# s<-summary(cv.modsel.rej)
# plot(cv.modsel.rej, names.arg=mods)

# Calculate posterior probabilities of each scenario (using mnlogistic)
modsel<-postpr(stats.obs, models, stat.sim2[,1:t], tol=0.001, method="mnlogistic")
modsel2<-postpr(stats.obs, models, stat.sim2[,1:t], tol=0.005, method="neuralnet")
modsel001<-postpr(stats.obs, models, stat.sim2[,1:t], tol=0.001, method="neuralnet")
summary(modsel)

# tol=0.001 produces infinite values for some models; trim to higher likelihood models
# and re-run.
stat.sim.trim<-subset(stat.sim2, models=="dn" | models=="up")
models.trim<-models$stat.sim.trim
modsel.trim2<-postpr(stats.obs, models.trim, stat.sim.trim[,1:t], tol=0.001, method="neuralnet")

# # Calculate posterior probabilities of each scenario (using rejection)
# modsel.rej<-postpr(stats.obs, models, stat.sim2[,1:36], tol=0.05, method="rejection")
# summary(modsel.rej)

# Goodness-of-fit
# Plot histogram of null distribution against best model and superimpose observed value
# Test goodness-of-fit
stat.sim.mod<-stat.sim2[stat.sim2$models=="am",1:t]
res.gfit.am<-gfit(target=stats.obs, sumstat=stat.sim.mod,
               statistic=median, tol = 0.001, nb.replicate=1000)
summary(res.gfit)
plot(res.gfit, main="cpup, nb=1000, tol=0.05, stat=median, Histogram under H0")

# Good to check that other models don't also provide a good fit, by running
# gfit with different models.

# Look at position in PC space; cprob=0.1 shows 90% envelope
gfitpca2(target=stats.obs, sumstat=stat.sim2[,1:t], 
         index=models, cprob=0.1, xlim=c(-10, 10), ylim=c(-10, 10))

# Infer parameter values --------------------------------------------------
# Select simulated parameters that correspond to the desired model
stat.sim.mod<-subset(stat.sim2, models=="dn")
par.sim.mod<-subset(par.sim2, models=="dn")

# e.g. infer migration rate
# can ABC estimate the migration rate?
cv.res<-cv4abc(param=par.sim.mod[,"MIGR"], sumstat=stat.sim.mod[,1:(ncol(stat.sim.mod)-1)], nval=100,
               tols=c(0.05, 0.1), method="rejection")
# Collins et al 2014, PRSB, 281, 20140097 use tolerance values of 0.001, 0.002, and 0.004
# Prost et al 2013 GCB, 19, 1854-1864 has really well-described methods. 

summary(cv.res)

par(mar=c(5, 3, 4, 0.5), cex=0.8)
plot(cv.res, caption="rejection")

res<-abc(target=stats.obs, param=par.sim.mod[,"MIGR"],
         sumstat=stat.sim.mod[,1:(ncol(stat.sim.up)-1)], tol=0.05, method="rejection")
#generate summaries of the posterior distributions incl. median and credible intervals
summary(res)
hist(res)

# Posterior predictive checks ---------------------------------------------

# 1. estimate posterior distributions of the parameters of each model
# 2. sample a set of N multivariate parameters from their posterior distribution
# 3. obtain sample from the distribution of the summary statistics a posteriori
# by simulating the data sets with the N sampled multivariate parameters using simulation
# software
