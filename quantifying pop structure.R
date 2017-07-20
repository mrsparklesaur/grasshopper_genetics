# load libraries
library(ggplot2)

# load data
data<-read.table("./sanguinipes_MAF02.txt", sep="\t", as.is=T, check.names=F)
data<-t(data)

# PCA ---------------------------------------------------------------------

#If missing values, impute with mean
# a<-is.na(data)
# if (any(a =TRUE)) {
#   data<-apply(data,1,function(x){ix <- which(is.na(x));x[ix] <- mean(x,na.rm=T); return(x)})
#   data<-t(data)
# }

genotypes.centered <- scale(data,center=T,scale=F)
genotypes.svd <- svd(genotypes.centered)

#####scree plot
plot(genotypes.svd$d^2/sum(genotypes.svd$d^2),xlab="PC",ylab="Fraction Variation Explained")

##### show %variation explained by first 6 PCs
pc1 <- paste("(",(100*(round((genotypes.svd$d^2/sum(genotypes.svd$d^2))[1],digits=4))),"%",")",sep="")
pc2 <- paste("(",(100*(round((genotypes.svd$d^2/sum(genotypes.svd$d^2))[2],digits=4))),"%",")",sep="")

##### Project markers onto the PCs
genotypes.PCbasis <- genotypes.centered %*% genotypes.svd$v

#####Create dataframe with metadata
metadata <- read.table("./grasshoppers/sanguinipes_sites.txt",sep="\t",comment.char="",as.is=T,check.names=F, header=T)
metadata2 <- metadata[match(rownames(data),metadata[,1]),]

##### plot full dataset
mainlab<-"M. boulderensis, 5501 markers, 66 individuals,\n MAF > 0.02, no missing data, impute"

# select colour pallette
library(RColorBrewer)

bg<-brewer.pal(5, "YlGnBu")

plot(genotypes.PCbasis[,3:4],col="black",bg=bg[factor(metadata2[,5])],
     xlab=paste("PC1",pc1,sep=" "),ylab=paste("PC2",pc2,sep=" "),
     pch=21,cex=1.5, main=mainlab, cex.main = 0.5)
legend("bottomleft",legend=c("RF","A1","B1","D1"),pt.bg=bg, pt.cex = 1.5, pch=21,cex=1,ncol=3)


datset<-cbind(data.frame("PC1" = genotypes.PCbasis[,1]), 
              data.frame("PC2" = genotypes.PCbasis[,2]), 
              data.frame("pop" = factor(metadata2[,5])))

a<-ggplot(datset, aes(PC1, PC2, fill=pop))
figA<-a + geom_point(shape = 21, size=3)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour="black", size=10))+
  xlab(paste("PC1 ",pc1, sep=""))+
  ylab(paste("PC2 ", pc2, sep=""))+
  scale_fill_manual(values = c(bg[1:5], bg[5]))

png("sanguinipes__PCA.png",
    type="cairo",
    width=100,
    height=80,
    units="mm",
    res=300)

figA

dev.off()
